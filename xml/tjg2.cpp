#include "FarmDb.hpp"
#include "FarmGeo.hpp"
#include "BoundarySwaths.hpp"

#include <pugixml.hpp>

#include <exception>
#include <stdexcept>
#include <format>
#include <iostream>
#include <iomanip>

namespace rng = std::ranges;

int main(int argc, const char *argv[]) {
  try {
    const auto path = (argc > 1) ? argv[1] : "TASKDATA.XML";

    static const char* RootName = "ISO11783_TaskData";

    auto doc = pugi::xml_document{};
    auto res = doc.load_file(path, pugi::parse_default | pugi::parse_ws_pcdata);
    if (!res) {
      std::cerr << std::format("XML parse error: {} (offset {})\n",
                               res.description(), res.offset);
      return 1;
    }

    auto root = doc.child(RootName);
    if (!root) {
      std::cerr << std::format("error: missing root <{}>\n", RootName);
      return 1;
    }

    auto db = farm_db::FarmDb{root};
    std::cout << db.customers.size() << " customers\n"
              << db.farms.size()     << " farms\n"
              << db.fields.size()    << " fields\n"
              << db.values.size()    << " values\n";

    db.swVendor = "Terry Golubiewski";
    db.swVersion = "0.1 (alpha)";

    int nextGuideId = 1;
    int nextSwathId = 1;
    for (auto& field: db.fields) {
      int partNum = 0;
      for (const auto& part: field.parts) {
        auto guideId   = std::format("GGP{:02}", nextGuideId++);
        auto guideName = field.name + std::format("_{:02}", ++partNum);
        field.guides.emplace_back(guideId, guideName);
        auto& guide = field.guides.back();
        auto geoSwaths =
            farm_db::BoundarySwaths(farm_db::Geo(part), 20.0 * farm_db::metre);
        int swathNum = 0;
        for (const auto& geoSwath: geoSwaths) {
          auto paths = std::vector<farm_db::Path>{};
          paths.reserve(geoSwath.size());
          for (const auto& geoPath: geoSwath)
            paths.emplace_back(farm_db::MakePath(geoPath));
          auto swathId = std::format("GPN{:02}", nextSwathId++);
          guide.swaths.emplace_back(swathId, paths);
          guide.swaths.back().name =
                                  guideName + std::format("_{:02}", ++swathNum);
        }
      }
    }

    auto doc2 = pugi::xml_document{};
    auto root2 = doc2.append_child(RootName);
    db.dump(root2);

    auto ok = doc2.save_file("out2.xml", "  ");
    if (!ok)
      throw std::runtime_error{"Error writing 'out2.xml'"};
    return 0;
  }
  catch (std::exception& x) {
    std::cerr << "Exception: " << x.what() << '\n';
  }
  catch (...) {
    std::cerr << "Unknown exception\n";
  }

  return 1;
} // main

