// taskdata_parse.cpp
//
// Build (example):
//   g++ -std=gnu++23 -O2 -Wall -Wextra -pedantic \
//       taskdata_parse.cpp -lpugixml -o taskdata_parse
//
// Run:
//   ./taskdata_parse TASKDATA.XML
//
// Lines wrapped to <= 80 cols.

#include "get_attr.hpp"

#include <pugixml.hpp>

#include <unordered_map>
#include <string>
#include <vector>
#include <string_view>
#include <utility>
#include <cstdio>
#include <cstring>

// ---------------------------------------------------------------------
// ISOXML constants (centralized so parse/write share the same strings).
namespace isoxml {

inline constexpr char Root[] = "ISO11783_TaskData";
inline constexpr char CTR[]  = "CTR";
inline constexpr char FRM[]  = "FRM";
inline constexpr char PFD[]  = "PFD";
inline constexpr char PLN[]  = "PLN";
inline constexpr char LSG[]  = "LSG";
inline constexpr char PNT[]  = "PNT";

namespace root_attr {
  inline constexpr char VersionMajor[]        = "VersionMajor";
  inline constexpr char VersionMinor[]        = "VersionMinor";
  inline constexpr char DataTransferOrigin[]  = "DataTransferOrigin";
  inline constexpr char MgmtSoftwareManufacturer[] =
                                              "ManagementSoftwareManufacturer";
  inline constexpr char MgmtSoftwareVersion[] = "ManagementSoftwareVersion";
} // root_attr

} // isoxml

// ---------------------------------------------------------------------
// Minimal in-memory shapes for a demo summary.
struct RootMeta {
  int versionMajor       = -1;  // required
  int versionMinor       = -1;  // required
  int dataTransferOrigin = -1;  // optional, -1 means "unset"
  std::string mgmtManufacturer; // optional
  std::string mgmtVersion;      // optional
  std::vector<std::pair<std::string,std::string>> unknown;
}; // RootMeta

struct Customer {
  std::string id;
  std::string name;
};

struct Farm {
  std::string id;
  std::string name;
  std::string ctrId;
};

struct Field {
  std::string id;
  std::string name;
  std::string ctrId;
  std::string frmId;
  // Minimal geometry stats for a quick sanity print.
  int polygonCount = 0;
  int ringCount    = 0; // LSG
  int pointCount   = 0; // PNT
};

struct Model {
  RootMeta root;
  std::unordered_map<std::string, Customer> customers;
  std::unordered_map<std::string, Farm>     farms;
  std::unordered_map<std::string, Field>    fields;
};

// ---------------------------------------------------------------------
// Diagnostics.
namespace {

void WarnUnknownAttr(const char *elem, const char *key)
  { std::fprintf(stderr, "warning: <%s> unknown @%s\n", elem, key); }
void WarnUnknownElem(const char *parent, const char *child) {
  std::fprintf(stderr, "warning: <%s> has unknown child <%s>\n",
               parent, child);
}
void ErrorMissing(const char *elem, const char *key)
  { std::fprintf(stderr, "error: <%s> missing @%s\n", elem, key); }
void ErrorRef(const char *elem, const char *msg)
  { std::fprintf(stderr, "error: <%s> %s\n", elem, msg); }

} // local

// ---------------------------------------------------------------------
// Parse root meta (pugixml only; no custom helpers).
RootMeta ParseRootMeta(const pugi::xml_node &root) {
  RootMeta m;
  for (const auto& a : root.attributes()) {
    if (!a)
      break;
    const auto k = a.name();
    if (!k || !k[0])
      break;
    if (     std::strcmp(k, isoxml::root_attr::VersionMajor) == 0)
      m.versionMajor       = tjg::get_attr<int>(a);
    else if (std::strcmp(k, isoxml::root_attr::VersionMinor) == 0)
      m.versionMinor       = tjg::get_attr<int>(a);
    else if (std::strcmp(k, isoxml::root_attr::DataTransferOrigin) == 0)
      m.dataTransferOrigin = tjg::get_attr<int>(a);
    else if (std::strcmp(k, isoxml::root_attr::MgmtSoftwareManufacturer) == 0)
      m.mgmtManufacturer   = tjg::get_attr<std::string>(a);
    else if (std::strcmp(k, isoxml::root_attr::MgmtSoftwareVersion) == 0)
      m.mgmtVersion       = tjg::get_attr<std::string>(a);
    else
      m.unknown.emplace_back(k, a.value());
  }
  if (m.versionMajor < 0 || m.versionMinor < 0)
    std::fprintf(stderr, "error: root missing VersionMajor/VersionMinor\n");
  return m;
} // ParseRootMeta

// ---------------------------------------------------------------------
// CTR
void ParseCtr(const pugi::xml_node &n, Model &out) {
  auto custId   = std::optional<std::string>{};
  auto custName = std::optional<std::string>{};
  for (const auto& a : n.attributes()) {
    if (!a)
      break;
    const auto k = a.name();
    if (!k || !k[0])
      break;
    if (     std::strcmp(k, "A") == 0)
      custId   = tjg::get_attr<std::string>(a);
    else if (std::strcmp(k, "B") == 0)
      custName = tjg::get_attr<std::string>(a);
    else
      WarnUnknownAttr("CTR", k);
  }

  bool ok = true;
  if (!custId  ) { ErrorMissing("CTR", "A"); ok = false; }
  if (!custName) { ErrorMissing("CTR", "B"); ok = false; }
  if (!ok) return;

  Customer ctr;
  ctr.id   = *custId;
  ctr.name = std::move(*custName);
  auto ins = out.customers.try_emplace(std::move(*custId), std::move(ctr));
  if (!ins.second)
    ErrorRef("CTR", "duplicate id");
} // ParseCtr

// ---------------------------------------------------------------------
// FRM
void ParseFrm(const pugi::xml_node &n, Model &out) {
  auto custId   = std::optional<std::string>();
  auto farmId   = std::optional<std::string>();
  auto farmName = std::optional<std::string>();

  for (const auto& a : n.attributes()) {
    if (!a)
      break;
    const auto k = a.name();
    if (!k || !k[0])
      break;
    if (     std::strcmp(k, "A") == 0)
      farmId   = tjg::get_attr<std::string>(a);
    else if (std::strcmp(k, "B") == 0)
      farmName = tjg::get_attr<std::string>(a);
    else if (std::strcmp(k, "I") == 0)
      custId   = tjg::get_attr<std::string>(a);
    else
      WarnUnknownAttr("FRM", k);
  }

  bool ok = true;
  if (!farmId  ) { ErrorMissing("FRM", "A"); ok = false; }
  if (!farmName) { ErrorMissing("FRM", "B"); ok = false; }
  if (!custId  ) { ErrorMissing("FRM", "I"); ok = false; }
  if (!ok) return;

  Farm frm;
  frm.id    = *farmId;
  frm.name  = std::move(*farmName);
  frm.ctrId = std::move(*custId);

  if (!out.customers.contains(frm.ctrId))
    ErrorRef("FRM", "references missing CTR (may appear later)");

  auto ins = out.farms.emplace(std::move(*farmId), std::move(frm));
  if (!ins.second)
    ErrorRef("FRM", "duplicate id");
} // ParseFrm

// ---------------------------------------------------------------------
// Minimal PNT counting for a field.
void ParsePnt(const pugi::xml_node &pnt, Field &F) {
  ++F.pointCount;
  for (const auto& a : pnt.attributes()) {
    if (!a)
      break;
    const auto k = a.name();
    if (!k || !k[0])
      break;
    if (   std::strcmp(k, "A") && std::strcmp(k, "B")
        && std::strcmp(k, "C") && std::strcmp(k, "D")) {
      WarnUnknownAttr(isoxml::PNT, k);
    }
  }
  for (const auto& c : pnt.children()) {
    if (!c)
      break;
    const auto k = c.name();
    if (!k || !k[0])
      break;
    WarnUnknownElem(isoxml::PNT, k);
  }
} // ParsePnt

// ---------------------------------------------------------------------
// Minimal LSG counting for a field.
void ParseLsg(const pugi::xml_node &lsg, Field &F) {
  ++F.ringCount;
  for (const auto& a : lsg.attributes()) {
    if (!a)
      break;
    const auto k = a.name();
    if (!k || !k[0])
      break;
    if (   std::strcmp(k, "A") && std::strcmp(k, "B")
        && std::strcmp(k, "C") && std::strcmp(k, "D")) {
      WarnUnknownAttr(isoxml::LSG, k);
    }
  }
  for (const auto& c : lsg.children()) {
    if (!c)
      break;
    const auto k = c.name();
    if (!k || !k[0])
      break;
    // std::printf("TJG '%s'\n", k);
    if (std::strcmp(k, isoxml::PNT) != 0) {
      WarnUnknownElem(isoxml::LSG, k);
      continue;
    }
    ParsePnt(c, F);
  }
} // ParseLsg

// ---------------------------------------------------------------------
// Minimal PLN parsing.
void ParsePln(const pugi::xml_node &pln, Field &F) {
  ++F.polygonCount;
  for (const auto& a : pln.attributes()) {
    if (!a)
      break;
    const auto k = a.name();
    if (!k || !k[0])
      break;
    if (std::strcmp(k, "A") && std::strcmp(k, "B"))
      WarnUnknownAttr(isoxml::PLN, k);
  }
  for (const auto& c : pln.children()) {
    if (!c)
      break;
    const auto k = c.name();
    if (!k || !k[0])
      break;
    if (std::strcmp(k, isoxml::LSG) != 0) {
      WarnUnknownElem(isoxml::PLN, k);
      continue;
    }
    ParseLsg(c, F);
  }
} // ParsePln

// ---------------------------------------------------------------------
// PFD
void ParsePfd(const pugi::xml_node &n, Model &out) {
  auto fieldId   = std::optional<std::string>{};
  auto fieldName = std::optional<std::string>{};
  auto ctrId     = std::optional<std::string>{};
  auto frmId     = std::optional<std::string>{};

  // Attributes (A,C,E,F)
  for (const auto& a : n.attributes()) {
    if (!a)
      break;
    const auto k = a.name();
    if (!k || !k[0])
      break;
    if (     std::strcmp(k, "A") == 0)
      fieldId   = tjg::get_attr<std::string>(a);
    else if (std::strcmp(k, "C") == 0)
      fieldName = tjg::get_attr<std::string>(a);
    else if (std::strcmp(k, "E") == 0)
      ctrId     = tjg::get_attr<std::string>(a);
    else if (std::strcmp(k, "F") == 0)
      frmId     = tjg::get_attr<std::string>(a);
    else
      WarnUnknownAttr("PFD", k);
  }

  bool ok = true;
  if (!fieldId  ) { ErrorMissing("PFD", "A"); ok = false; }
  if (!fieldName) { ErrorMissing("PFD", "B"); ok = false; }
  if (!ctrId    ) { ErrorMissing("PFD", "E"); ok = false; }
  if (!frmId    ) { ErrorMissing("PFD", "F"); ok = false; }
  if (!ok) return;

  Field f;
  f.id    = std::move(*fieldId);
  f.name  = std::move(*fieldName);
  f.ctrId = std::move(*ctrId);
  f.frmId = std::move(*frmId);

  // Children in doc order. We only count PLN/LSG/PNT; warn for others.
  for (const auto& c : n.children()) {
    if (!c)
      break;
    const auto k = c.name();
    if (!k || !k[0])
      break;
    if (std::strcmp(k, isoxml::PLN) != 0) {
      WarnUnknownElem("PFD", k);
      continue;
    }
    ParsePln(c, f);
  }

  if (!out.customers.contains(f.ctrId))
    ErrorRef("PFD", "references missing CTR (may appear earlier/later)");
  if (!out.farms.contains(f.frmId))
    ErrorRef("PFD", "references missing FRM (may appear earlier/later)");

  auto id = f.id;
  auto ins = out.fields.emplace(id, std::move(f));
  if (!ins.second)
    ErrorRef("PFD", "duplicate id");
} // ParsePfd

// ---------------------------------------------------------------------
// Pretty summary.
void PrintSummary(const Model &m) {
  std::printf("Root: Version %d.%d, DTO=%d\n",
              m.root.versionMajor, m.root.versionMinor,
              m.root.dataTransferOrigin);
  if (!m.root.mgmtManufacturer.empty() || !m.root.mgmtVersion.empty()) {
    std::printf("Mgmt: %s %s\n",
      m.root.mgmtManufacturer.c_str(),
      m.root.mgmtVersion.c_str());
  }
  if (!m.root.unknown.empty()) {
    std::printf("Root unknown attrs:\n");
    for (const auto &kv : m.root.unknown)
      std::printf("  %s=\"%s\"\n", kv.first.c_str(), kv.second.c_str());
  }

  std::printf("Counts: CTR=%zu, FRM=%zu, PFD=%zu\n",
              m.customers.size(), m.farms.size(), m.fields.size());

  int poly = 0, rings = 0, pts = 0;
  for (const auto &kv : m.fields) {
    poly  += kv.second.polygonCount;
    rings += kv.second.ringCount;
    pts   += kv.second.pointCount;
  }
  std::printf("Geometry: PLN=%d, LSG=%d, PNT=%d\n", poly, rings, pts);
} // PrintSummary

// ---------------------------------------------------------------------
// Main: load, basic validation, ordered traversal with dispatch.
int main(int argc, const char *argv[]) {
  const auto path = (argc > 1) ? argv[1] : "TASKDATA.XML";

  pugi::xml_document doc;
  auto res = doc.load_file(path, pugi::parse_default | pugi::parse_ws_pcdata);
  if (!res) {
    std::fprintf(stderr, "XML parse error: %s (offset %zu)\n",
                 res.description(), static_cast<size_t>(res.offset));
    return 1;
  }

  auto root = doc.child(isoxml::Root);
  if (!root) {
    std::fprintf(stderr, "error: missing root <%s>\n", isoxml::Root);
    return 1;
  }

  Model model;
  model.root = ParseRootMeta(root);

  // Walk children in file order. Dispatch on name via strcmp.
  for (const auto& n : root.children()) {
    if (!n)
      break;
    const auto k = n.name();
    if (!k || !k[0])
      break;
    if (     std::strcmp(k, isoxml::CTR) == 0)
      ParseCtr(n, model);
    else if (std::strcmp(k, isoxml::FRM) == 0)
      ParseFrm(n, model);
    else if (std::strcmp(k, isoxml::PFD) == 0)
      ParsePfd(n, model);
    else
      WarnUnknownElem(isoxml::Root, k);
  }

  PrintSummary(model);
  return 0;
} // main

