#include "FarmDb.hpp"

#include <pugixml.hpp>

#include <exception>
#include <stdexcept>
#include <iostream>
#include <iomanip>

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

