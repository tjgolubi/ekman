/// @file
/// @brief Read ISO11783 TaskData.xml (or similar) and extract CTR/FRM and
///        PFD/PLN/LSG/PNT data. Writes per-LSG .xy when --dump is used.
/// Usage:
///   field_xml <taskdata.xml> [--dump] [--outdir DIR]
///
/// Conventions:
/// - C++23, almost-always-auto, exceptions, fs::path, fstreams, anon ns.
/// - PNT attributes: C=latitude, D=longitude (WGS84 degrees).
/// - Filenames sanitized with isalnum + [_-.].
/// - No schema validation; fast DOM parsing via pugixml.

#include <gsl-lite/gsl-lite.hpp>

#include <pugixml.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <cctype>

namespace {

namespace fs = std::filesystem;

// ---------- Data types ----------

struct Pt {
  double lon = 0.0;  // x
  double lat = 0.0;  // y
};

struct Lsg {
  std::string pfd_id;   // owning field id
  std::string pln_id;   // plan id
  std::string lsg_id;   // this group id
  std::vector<Pt> pts;  // (lon,lat) degrees
};

struct Customer {
  std::string id;    // CTR @A
  std::string name;  // CTR @B
};

struct Farm {
  std::string id;            // FRM @A
  std::string name;          // FRM @B
  std::string customer_id; // FRM @I -> CTR @A
};

struct FieldMeta {
  // From PFD
  std::string id;           // PFD @A
  std::string name;         // PFD @C (optional display)
  std::string customer;   // PFD @E -> CTR @A
  std::string farm;       // PFD @F -> FRM @A
};

// ---------- Helpers ----------

auto Sanitize(std::string_view s) -> std::string {
  auto out = std::string{};
  out.reserve(s.size());
  for (auto ch : s) {
    const auto c = static_cast<unsigned char>(ch);
    out.push_back(std::isalnum(c) || c == '_' || c == '-' || c == '.'
                    ? static_cast<char>(c)
                    : '_');
  }
  if (out.empty()) out = "na";
  return out;
}

auto RequireAttr(const pugi::xml_node& n, const char* key) -> std::string {
  auto a = n.attribute(key);
  if (a) {
    auto s = std::string{a.value()};
    return s;
  }
  throw std::runtime_error{
      std::string{"missing attribute \""} + key + "\" on <" + n.name() + ">"};
}

// ---------- Parse ----------

struct Parsed {
  // Primary entities
  std::unordered_map<std::string, Customer> customers; // by CTR @A
  std::unordered_map<std::string, Farm>     farms;     // by FRM @A
  std::unordered_map<std::string, FieldMeta>  fields;      // by PFD @A

  // Geometries
  std::vector<Lsg> lsgs; // all PFD/PLN/LSG rings (lon,lat)
}; // Parsed

auto LoadTaskData(const fs::path& xml_path) -> Parsed {
  auto doc = pugi::xml_document{};

  // Keep a stable string; pugixml's load_file expects c-string
  const auto path_str = xml_path.string();
  const auto res = doc.load_file(path_str.c_str());
  if (!res) {
    auto msg = std::string{"XML parse error for '"} + xml_path.string()
             + "': " + res.description();
    throw std::runtime_error(msg);
  }

  const auto root = doc.child("ISO11783_TaskData");
  if (!root)
    throw std::runtime_error("not an <ISO11783_TaskData> document");

  auto out = Parsed{};

  // --- Customers (CTR) ---
  for (auto ctr : root.children("CTR")) {
    auto id   = RequireAttr(ctr, "A");
    auto name = Sanitize(RequireAttr(ctr, "B"));
    out.customers.emplace(id, Customer{id, std::move(name)});
  }

  // --- Farms (FRM) ---
  for (auto frm : root.children("FRM")) {
    auto id  = RequireAttr(frm, "A");
    auto nam = Sanitize(RequireAttr(frm, "B"));
    auto cid = RequireAttr(frm, "I"); // customer id
    out.farms.emplace(id, Farm{id, std::move(nam), std::move(cid)});
  }

  // --- Fields (PFD) + Plans + LineStringGroups ---
  for (auto pfd : root.children("PFD")) {
    // Field meta (ids may be empty in some exports; treat as "")
    auto pfd_id   = std::string{pfd.attribute("A").as_string("")};
    auto pfd_name = std::string{pfd.attribute("C").as_string("")};
    auto ctr_id   = std::string{pfd.attribute("E").as_string("")}; // customer ref
    auto frm_id   = std::string{pfd.attribute("F").as_string("")}; // farm ref

    if (!pfd_id.empty()) {
      out.fields.emplace(pfd_id,
                         FieldMeta{pfd_id, pfd_name, ctr_id, frm_id});
    }

    for (auto pln : pfd.children("PLN")) {
      auto pln_id = std::string{pln.attribute("A").as_string("")};

      for (auto lsg : pln.children("LSG")) {
        auto ring = Lsg{};
        ring.pfd_id = pfd_id;
        ring.pln_id = pln_id;
        ring.lsg_id = std::string{lsg.attribute("A").as_string("")};

        for (auto pnt : lsg.children("PNT")) {
          const auto lat = pnt.attribute("C").as_double(
              std::numeric_limits<double>::quiet_NaN());
          const auto lon = pnt.attribute("D").as_double(
              std::numeric_limits<double>::quiet_NaN());
          if (!std::isfinite(lat) || !std::isfinite(lon))
            throw std::runtime_error("invalid <PNT> coordinates");
          ring.pts.push_back(Pt{lon, lat}); // store lon=x, lat=y
        }

        if (!ring.pts.empty())
          out.lsgs.push_back(std::move(ring));
      }
    }
  }

  // Optional: basic referential checks (don’t hard-fail; just warn)
  for (const auto& [fid, f] : out.fields) {
    if (!f.customer.empty() && !out.customers.contains(f.customer)) {
      std::cerr << "warning: PFD " << fid
                << " references missing CTR " << f.customer << '\n';
    }
    if (!f.farm.empty() && !out.farms.contains(f.farm)) {
      std::cerr << "warning: PFD " << fid
                << " references missing FRM " << f.farm << '\n';
    }
  }

  // Require *something* parsed
  if (out.customers.empty() && out.farms.empty()
      && out.fields.empty() && out.lsgs.empty())
    throw std::runtime_error("no recognizable ISO11783 data found");

  return out;
} // LoadTaskData

// ---------- Output ----------

void PrintSummary(const Parsed& data) {
  std::cout << "Customers: " << data.customers.size() << '\n';
  for (const auto& [id, c] : data.customers)
    std::cout << "  CTR " << id << " : " << c.name << '\n';

  std::cout << "Farms: " << data.farms.size() << '\n';
  for (const auto& [id, f] : data.farms)
    std::cout << "  FRM " << id << " : " << f.name
              << " (CTR=" << f.customer_id << ")\n";

  std::cout << "Fields: " << data.fields.size() << '\n';
  for (const auto& [fid, f] : data.fields) {
    const auto& cname = (f.customer.empty() || !data.customers.contains(f.customer))
                        ? std::string_view{"?"}
                        : std::string_view{data.customers.at(f.customer).name};
    const auto& fname = (f.farm.empty() || !data.farms.contains(f.farm))
                        ? std::string_view{"?"}
                        : std::string_view{data.farms.at(f.farm).name};
    std::cout << "  PFD " << fid
              << " : \"" << f.name << "\""
              << " (CTR=" << f.customer << ":" << cname
              << ", FRM=" << f.farm << ":" << fname << ")\n";
  }

  std::cout << "LineStringGroups (LSG): " << std::ssize(data.lsgs) << '\n';
  for (const auto& r : data.lsgs) {
    std::cout << "  PFD=" << r.pfd_id
              << " PLN=" << r.pln_id
              << " LSG=" << r.lsg_id
              << " pts=" << std::ssize(r.pts) << '\n';
  }
} // PrintSummary

void DumpRings(const std::vector<Lsg>& rings,
               const fs::path& outdir = fs::path{})
{
  for (const auto& r : rings) {
    const auto base = "pfd-"  + Sanitize(r.pfd_id)
                    + "_pln-" + Sanitize(r.pln_id)
                    + "_lsg-" + Sanitize(r.lsg_id) + ".xy";
    const auto fn = outdir / base;

    auto ofs = std::ofstream{fn, std::ios::binary};
    if (!ofs)
      throw std::runtime_error("cannot open for writing: " + fn.string());

    ofs.setf(std::ios::fixed, std::ios::floatfield);
    ofs.precision(12);
    for (const auto& p : r.pts)
      ofs << p.lon << ' ' << p.lat << '\n';

    std::cout << "wrote " << fn << " (" << std::ssize(r.pts) << " points)\n";
  }
} // DumpRings

} // local

// ---------- main ----------

auto main(int argc, char** argv) -> int {
  if (argc < 2) {
    std::cerr << "usage: field_xml <taskdata.xml> [--dump] [--outdir DIR]\n";
    return EXIT_FAILURE;
  }
  try {
    const auto xml = std::filesystem::path{argv[1]};

    auto do_dump = false;
    auto outdir  = std::filesystem::path{};

    for (int i = 2; i < argc; ++i) {
      const auto arg = std::string_view{argv[i]};
      if (arg == "--dump") {
        do_dump = true;
      } else if (arg == "--outdir" && i + 1 < argc) {
        outdir = std::filesystem::path{argv[++i]};
      }
    }

    const auto parsed = LoadTaskData(xml);
    PrintSummary(parsed);
    if (do_dump) DumpRings(parsed.lsgs, outdir);

    return EXIT_SUCCESS;
  }
  catch (const std::exception& e) {
    std::cerr << "std::exception: " << e.what() << '\n';
  }
  catch (...) {
    std::cerr << "unknown exception\n";
  }
  return EXIT_FAILURE;
} // main

