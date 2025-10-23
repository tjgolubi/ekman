// taskdata_parse.cpp
//
// Build (example):
//   g++ -std=gnu++23 -O2 -Wall -Wextra -pedantic
//       taskdata_parse.cpp -lpugixml -o taskdata_parse
//
// Run:
//   ./taskdata_parse TASKDATA.XML
//
// Lines wrapped to <= 80 cols.

#include "get_attr.hpp"

#include <pugixml.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/is_valid.hpp>
#include <boost/geometry/core/cs.hpp>

#include <unordered_map>
#include <regex>
#include <string>
#include <vector>
#include <string_view>
#include <functional>
#include <format>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <exception>
#include <utility>
#include <cstdint>

namespace ggl = boost::geometry;
using GeoPoint = ggl::model::point<double, 2, ggl::cs::geographic<ggl::degree>>;
using GeoLineString = ggl::model::linestring<GeoPoint>;
using GeoPolyLine   = ggl::model::multi_linestring<GeoLineString>;
using GeoRing       = ggl::model::ring<GeoPoint, true >;
using GeoHole       = ggl::model::ring<GeoPoint, false>;
using GeoPolygon    = ggl::model::polygon<GeoPoint>;

// ---------------------------------------------------------------------
// ISOXML constants (centralized so parse/write share the same strings).
namespace isoxml {

constexpr char Root[] = "ISO11783_TaskData";

const auto CtrRe = std::regex{"CTR-?[0-9]+"};
const auto FrmRe = std::regex{"FRM-?[0-9]+"};
const auto PfdRe = std::regex{"PFD-?[0-9]+"};
const auto GpnRe = std::regex{"GPN-?[0-9]+"};
const auto GgpRe = std::regex{"GGP-?[0-9]+"};

namespace root_attr {
  constexpr char VersionMajor[]        = "VersionMajor";
  constexpr char VersionMinor[]        = "VersionMinor";
  constexpr char DataTransferOrigin[]  = "DataTransferOrigin";
  constexpr char MgmtSoftwareManufacturer[] = "ManagementSoftwareManufacturer";
  constexpr char MgmtSoftwareVersion[] = "ManagementSoftwareVersion";
} // root_attr

} // isoxml

using XmlNode = pugi::xml_node;
using XmlAttr = pugi::xml_attribute;

[[noreturn]] void InvalidNode(const XmlNode& xml, std::string what) {
  auto k = tjg::name(xml);
  what.reserve(what.size() + 8 + k.size());
  what += " on <";
  what += k;
  what += ">";
  throw std::runtime_error{what};
} // InvalidNode

[[noreturn]] void InvalidAttr(const XmlNode& xml, const char* key,
                              std::string what="Invalid attribute")
{
  auto k = std::string_view{key};
  what.reserve(what.size() + k.size() + 40);
  what += " \"";
  what += k;
  what += "\" ";
  auto a = xml.attribute(k);
  if (a) {
    what += "= ";
    what += a.as_string();
  }
  else {
    what += "is missing";
  }
  InvalidNode(xml, what);
} // InvalidAttr

template<typename T>
std::optional<T> GetAttr(const XmlNode& x, const char* key) {
  const auto& a = x.attribute(key);
  return (a) ? tjg::try_get_attr<T>(a) : std::nullopt;
} // GetAttr

template<typename T=std::string>
T RequireAttr(const pugi::xml_node& x, const char* key) {
  const auto& a = x.attribute(key);
  if (!a)
    InvalidAttr(x, key);
  auto v = GetAttr<T>(x, key);
  if (v)
    return *v;
  if constexpr (std::is_same_v<T, std::string>) {
    auto s = a.as_string();
    if (s[0] != '\0')
      return s;
  }
  InvalidAttr(x, key);
} // RequireAttr

struct Customer {    // CTR
  std::string id;    // A
  std::string name;  // B
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  Customer() = default;
  Customer(std::string_view id_, std::string_view name_)
    : id{id_}, name{name_}
  {
    if (!std::regex_match(id, isoxml::CtrRe))
      throw std::invalid_argument{"Customer: invalid id " + id};
  }
  explicit Customer(const XmlNode& node);
  void dump(XmlNode& node) const;
}; // Customer

Customer::Customer(const XmlNode& node)
  : id  {RequireAttr<std::string>(node, "A")}
  , name{RequireAttr<std::string>(node, "B")}
{
  for (const auto& a: node.attributes()) {
    auto k = tjg::name(a);
    if (k == "A" || k == "B")
      continue;
    otherAttrs.emplace_back(a.name(), a.value());
  }
} // Customer::ctor

void Customer::dump(XmlNode& node) const {
  node.set_name("CTR");
  node.append_attribute("A") = id;
  node.append_attribute("B") = name;
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
} // Customer::dump

struct Farm {           // FRM
  std::string id;       // A
  std::string name;     // B
  std::string ctrId;    // I optional
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  Farm() = default;
  Farm(std::string_view id_, std::string_view name_)
    : id{id_}, name{name_}
  {
    if (!std::regex_match(id, isoxml::FrmRe))
      throw std::invalid_argument{"Farm: invalid id: " + id};
  }
  explicit Farm(const XmlNode& node);
  void dump(XmlNode& node) const;
}; // Farm

Farm::Farm(const XmlNode& node)
  : id   {RequireAttr<std::string>(node, "A")}
  , name {RequireAttr<std::string>(node, "B")}
{
  for (const auto& a: node.attributes()) {
    auto k = tjg::name(a);
    if (k == "A" || k == "B")
      continue;
    if (k == "I") ctrId = tjg::get_attr<std::string>(a);
    else otherAttrs.emplace_back(a.name(), a.value());
  }
} // Farm::ctor

void Farm::dump(XmlNode& node) const {
  node.set_name("FRM");
  node.append_attribute("A") = id;
  node.append_attribute("B") = name;
  node.append_attribute("I") = ctrId;
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
} // Farm::dump

using Angle = double;
using LonLat = GeoPoint;
Angle Latitude (const LonLat& ll) { return ggl::get<1>(ll); }
Angle Longitude(const LonLat& ll) { return ggl::get<0>(ll); }
LonLat LatLon(double lat, double lon) {
  if ( lat <  -90.0 || lat >  +90.0
    || lon < -180.0 || lon > +180.0)
  {
    auto msg = std::format("LatLong: invalid ({}, {})", lat, lon);
    throw std::range_error{msg};
  }
  return LonLat{lon, lat};
} // LatLon::ctor

struct Point {
  enum class Type {
    Flag=1, Other, Access, Storage, Obstacle, GuideA, GuideB,
    GuideCenter, GuidePoint, Field, Base
  };
  using Types = tjg::EnumList<Type, Type::Flag, Type::Other, Type::Access,
      Type::Storage, Type::Obstacle, Type::GuideA, Type::GuideB,
      Type::GuideCenter, Type::GuidePoint, Type::Field, Type::Base>;
  Type type;
  LonLat point;
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  using TypeValidator = std::function<bool(Type)>;
  static constexpr bool AnyType(Type) { return true; }
  GeoPoint geo() const { return point; }
  Angle latitude()  const { return Latitude(point);  }
  Angle longitude() const { return Longitude(point); }
  Point() = default;
  Point(Angle lat, Angle lon, Type type_ = Type::Other)
    : type{type_}, point{lon, lat} { }
  explicit Point(const GeoPoint& geo_pt, Type type_ = Type::Other)
    : type{type_}, point{geo_pt} { }
  explicit Point(const XmlNode& x, TypeValidator validator=AnyType);
  void dump(XmlNode& x) const;
}; // Point

const char* Name(Point::Type x) noexcept {
  switch (x) {
    case Point::Type::Flag:        return "Flag";
    case Point::Type::Other:       return "Other";
    case Point::Type::Access:      return "Access";
    case Point::Type::Storage:     return "Storage";
    case Point::Type::Obstacle:    return "Obstacle";
    case Point::Type::GuideA:      return "GuideA";
    case Point::Type::GuideB:      return "GuideB";
    case Point::Type::GuideCenter: return "GuideCenter";
    case Point::Type::GuidePoint:  return "GuidePoint";
    case Point::Type::Field:       return "Field";
    case Point::Type::Base:        return "Base";
    default: return nullptr;
  }
} // Name(Type)

namespace tjg {
template<> struct EnumValues<Point::Type>: Point::Types { };
} // tjg

Point::Point(const XmlNode& x, TypeValidator validator)
  : type{ RequireAttr<Type >(x, "A")}
  , point{LatLon(RequireAttr<Angle>(x, "C"),
                 RequireAttr<Angle>(x, "D"))}
{
  for (const auto& a: x.attributes()) {
    auto k = tjg::name(a);
    if (k == "A" || k == "C" || k == "D")
      continue;
    otherAttrs.emplace_back(k, a.value());
  }
  if (!validator(type)) {
    auto msg = std::string{"Point: invalid type in this context: "}
             + Name(type);
    throw std::range_error{msg};
  }
} // ctor

void Point::dump(XmlNode& x) const {
  x.set_name("PNT");
  x.append_attribute("A") = static_cast<int>(type);
  x.append_attribute("C") = Latitude(point);
  x.append_attribute("D") = Longitude(point);
  for (const auto& [k, v]: otherAttrs)
    x.append_attribute(k) = v;
} // Point::dump

struct LineString { // LSG
  enum class Type {
    Exterior=1,
    Interior,
    TramLine,
    Sampling,
    Guidance,
    Drainage,
    Fence,
    Flag,
    Obstacle
  };
  using Types = tjg::EnumList<Type, Type::Exterior, Type::Interior,
      Type::TramLine, Type::Sampling, Type::Guidance, Type::Drainage,
      Type::Fence, Type::Flag, Type::Obstacle>;
  Type type;
  std::vector<Point> points;
  std::vector<std::pair<std::string, std::string>> otherAttrs;

  bool empty() const { return points.empty(); }
  std::size_t size() const { return points.size(); }

  GeoLineString geo() const {
    auto out = GeoLineString{};
    out.reserve(points.size());
    for (const auto& p: points)
      out.emplace_back(p.geo());
    return out;
  } // geo

  GeoRing ring() const {
    auto ls = geo();
    auto out = GeoRing{ls.begin(), ls.end()};
    ggl::correct(out);
    auto msg = std::string{};
    if (!ggl::is_valid(out, msg))
      throw std::runtime_error{"LineString::ring: not a ring: " + msg};
    return out;
  } // ring

  GeoHole hole() const {
    auto ls = geo();
    auto out = GeoHole{ls.begin(), ls.end()};
    ggl::correct(out);
    auto msg = std::string{};
    if (!ggl::is_valid(out, msg))
      throw std::runtime_error{"LineString::hole: not a hole: " + msg};
    return out;
  } // hole

  LineString() = default;
  LineString(Type type_, Point::Type ptType, const std::vector<LonLat>& pts)
    : type{type_}
  {
    points.reserve(pts.size());
    for (const auto& p: pts)
      points.emplace_back(p, ptType);
  }

  explicit LineString(const XmlNode& x,
                      Point::TypeValidator point_validator=Point::AnyType);
  void dump(XmlNode& node) const;
}; // LineString

const char* Name(LineString::Type x) noexcept {
  switch (x) {
    case LineString::Type::Exterior:  return "Exterior";
    case LineString::Type::Interior:  return "Interior";
    case LineString::Type::TramLine:  return "TramLine";
    case LineString::Type::Sampling:  return "Sampling";
    case LineString::Type::Guidance:  return "Guidance";
    case LineString::Type::Drainage:  return "Drainage";
    case LineString::Type::Fence:     return "Fence";
    case LineString::Type::Flag:      return "Flag";
    case LineString::Type::Obstacle:  return "Obstacle";
    default: return nullptr;
  }
} // Name(LineString::Type)

namespace tjg {
template<> struct EnumValues<LineString::Type>: LineString::Types { };
} // tjg

LineString::LineString(const XmlNode& x, Point::TypeValidator point_validator)
  : type{RequireAttr<Type>(x, "A")}
{
  for (const auto& a: x.attributes()) {
    auto k = tjg::name(a);
    if (k == "A")
      continue;
    otherAttrs.emplace_back(k, a.value());
  }
  for (const auto& c: x.children()) {
    if (c.type() != pugi::node_element)
      continue;
    auto k = tjg::name(c);
    if (k == "PNT") {
      points.emplace_back(c, point_validator);
      continue;
    }
    std::cerr << "LineString: element ignored: " <<  k << '\n';
  }
} // LineString ctor

void LineString::dump(XmlNode& node) const {
  node.set_name("LSG");
  node.append_attribute("A") = static_cast<int>(type);
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
  for (const auto& p: points) {
    auto c = node.append_child("PNT");
    p.dump(c);
  }
} // LineString::dump

struct Polygon { // PLN
  enum class Type {
    Boundary=1, Treatment, Water, Building, Road, Obstacle, Flag, Other, Field,
    Headland, Buffer, Windbreak
  }; // Type

  using Types = tjg::EnumList<Type, Type::Boundary, Type::Treatment,
      Type::Water, Type::Building, Type::Road, Type::Obstacle, Type::Flag,
      Type::Other, Type::Field, Type::Headland, Type::Buffer, Type::Windbreak>;

  Type type;
  LineString outer;
  std::vector<LineString> inners;
  std::vector<std::pair<std::string, std::string>> otherAttrs;

  GeoPolygon geo() const {
    auto out = GeoPolygon{};
    out.outer() = outer.ring();
    for (const auto& p: inners)
      out.inners().emplace_back(p.ring());
    ggl::correct(out);
    auto msg = std::string{};
    if (!ggl::is_valid(out, msg))
      throw std::runtime_error{"Polygon::geo: invalid polygon: " + msg};
    return out;
  } // geo

  Polygon() = default;
  explicit Polygon(const GeoPolygon& poly, Type type_ = Type::Boundary);
  explicit Polygon(const XmlNode& x);
  void dump(XmlNode& node) const;
}; // Polygon

const char* Name(Polygon::Type x) noexcept {
  switch (x) {
    case Polygon::Type::Boundary:  return "Boundary";
    case Polygon::Type::Treatment: return "Treatment";
    case Polygon::Type::Water:     return "Water";
    case Polygon::Type::Building:  return "Building";
    case Polygon::Type::Road:      return "Road";
    case Polygon::Type::Obstacle:  return "Obstacle";
    case Polygon::Type::Flag:      return "Flag";
    case Polygon::Type::Other:     return "Other";
    case Polygon::Type::Field:     return "Field";
    case Polygon::Type::Headland:  return "Headland";
    case Polygon::Type::Buffer:    return "Buffer";
    case Polygon::Type::Windbreak: return "Windbreak";
    default: return nullptr;
  }
} // Name(Polygon::Type)

namespace tjg {
template<> struct EnumValues<Polygon::Type>: Polygon::Types { };
} // tjg

Polygon::Polygon(const GeoPolygon& poly, Type type_)
  : type{type_}
  , outer{LineString::Type::Exterior, Point::Type::Field, poly.outer()}
  , inners{}
{
  inners.reserve(poly.inners().size());
  for (const auto& r: poly.inners())
    inners.emplace_back(LineString::Type::Interior, Point::Type::Field, r);
}

Polygon::Polygon(const XmlNode& x)
  : type{RequireAttr<int>(x, "A")}
{
  auto point_validator = [](Point::Type t) -> bool {
    return (t == Point::Type::Field);
  };
  for (const auto& lsg: x.children("LSG")) {
    auto ring = LineString{lsg, point_validator};
    switch (ring.type) {
      case LineString::Type::Exterior:
        if (!outer.empty())
          throw std::runtime_error{"Polygon: multiple exterior rings"};
        outer = std::move(ring);
        break;
      case LineString::Type::Interior:
        inners.emplace_back(std::move(ring));
        break;
      default: {
        auto msg = std::string{"Polygon: unexpected LineString type: "}
                 + Name(ring.type);
        throw std::runtime_error{msg};
      }
    }
  }
  for (const auto& c: x.children()) {
    if (c.type() != pugi::node_element)
      continue;
    auto k = tjg::name(c);
    if (k == "LSG")
      continue;
    std::cerr << "Polygon: element ignored: " << k << '\n';
  }
  if (outer.empty())
    throw std::runtime_error{"Polygon: missing exterior ring"};
  if (std::ssize(outer) < 4)
    throw std::runtime_error{"Polygon: exterior ring too small"};
  for (const auto& r: inners) {
    if (std::ssize(r) < 4)
      throw std::runtime_error{"Polygon: inter ring too small"};
  }
} // Polygon ctor

void Polygon::dump(XmlNode& node) const {
  node.set_name("PLN");
  node.append_attribute("A") = static_cast<int>(type);
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
  auto c = node.append_child("LSG");
  outer.dump(c);
  for (const auto& lsg: inners) {
    c = (node.append_child("LSG"));
    lsg.dump(c);
  }
} // Polygon::dump

struct Swath { // GPN
  enum class Type { AB = 1, APlus, Curve, Pivot, Spiral };
  using TypeList = tjg::EnumList<Type, Type::AB, Type::APlus, Type::Curve,
                                 Type::Pivot, Type::Spiral>;
  enum class Option { CW=1, CCW, Full };
  using OptionList =
                  tjg::EnumList<Option, Option::CW, Option::CCW, Option::Full>;
  enum class Direction { Both=1, Left, Right, None };
  using DirectionList = tjg::EnumList<Direction, Direction::Both,
                            Direction::Left, Direction::Right, Direction::None>;
  enum class Extension { Both=1, First, Last, None };
  using ExtensionList = tjg::EnumList<Extension, Extension::Both,
                            Extension::First, Extension::Last, Extension::None>;
  enum class Method {
    NoGps=0, GNSS, DGNSS, PreciseGNSS, RtkInt, RtkFloat, DR, Manual, Sim,
    PC=16, Other
  };
  using MethodList = tjg::EnumList<Method, Method::NoGps, Method::GNSS,
      Method::DGNSS, Method::PreciseGNSS, Method::RtkInt, Method::RtkFloat,
      Method::DR, Method::Manual, Method::Sim, Method::PC, Method::Other>;
  std::string id;
  std::string name;       // optional
  Type type;
  std::optional<Option>    option;
  std::optional<Direction> direction;
  std::optional<Extension> extension;
  std::optional<Angle>     heading;
  std::optional<unsigned>  radius;
  std::optional<Method>    method;
  std::optional<double>    horizontalAccuracy;
  std::optional<double>    verticalAccuracy;
  std::string baseId;
  std::string srid;
  std::optional<unsigned>  leftRemaining;
  std::optional<unsigned>  rightRemaining;
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  std::vector<LineString> paths;
  std::vector<Polygon> boundaries;
  void push(const GeoLineString& path) {
    if (path.empty())
      return;
    auto ls = LineString{LineString::Type::Guidance, Point::Type::GuidePoint,
                         path};
    if (ls.size() > 1) {
      ls.points.front().type = Point::Type::GuideA;
      ls.points.back() .type = Point::Type::GuideB;
    }
    paths.emplace_back(std::move(ls));
  }
  Swath() = default;
  explicit Swath(std::string_view id_, Type type_ = Type::Curve)
    : id{id_}, type{type_}
  {
    if (!std::regex_match(id, isoxml::GpnRe))
      throw std::invalid_argument{"Swath: invalid id: " + id};
  }

  Swath(std::string_view id_, const GeoPolyLine& lines, Type type_=Type::Curve)
    : Swath{id_, type_}
  {
    for (const auto& l: lines)
      push(l);
  }

  Swath(std::string_view id_, const GeoLineString& ls, Type type_=Type::Curve)
    : Swath{id_, type_}
    { push(ls); }
  explicit Swath(const XmlNode& node);
  void dump(XmlNode& node) const;
}; // Swath

const char* Name(Swath::Type x) noexcept {
  switch (x) {
    case Swath::Type::AB:     return "AB";
    case Swath::Type::APlus:  return "APlus";
    case Swath::Type::Curve:  return "Curve";
    case Swath::Type::Pivot:  return "Pivot";
    case Swath::Type::Spiral: return "Spiral";
    default: return nullptr;
  }
} // Name(Swath::Type)

const char* Name(Swath::Option x) noexcept {
  switch (x) {
    case Swath::Option::CW:   return "CW";
    case Swath::Option::CCW:  return "CCW";
    case Swath::Option::Full: return "Full";
    default: return nullptr;
  }
} // Name(Swath::Option)

const char* Name(Swath::Direction x) noexcept {
  switch (x) {
    case Swath::Direction::Both:  return "Both";
    case Swath::Direction::Left:  return "Left";
    case Swath::Direction::Right: return "Right";
    case Swath::Direction::None:  return "None";
    default: return nullptr;
  }
} // Name(Swath::Direction)

const char* Name(Swath::Extension x) noexcept {
  switch (x) {
    case Swath::Extension::Both:  return "Both";
    case Swath::Extension::First: return "First";
    case Swath::Extension::Last:  return "Last";
    case Swath::Extension::None:  return "None";
    default: return nullptr;
  }
} // Name(Swath::Extension)

const char* Name(Swath::Method x) noexcept {
  switch (x) {
    case Swath::Method::NoGps:       return "NoGps";
    case Swath::Method::GNSS:        return "GNSS";
    case Swath::Method::DGNSS:       return "DGNSS";
    case Swath::Method::PreciseGNSS: return "PreciseGNSS";
    case Swath::Method::RtkInt:      return "RtkInt";
    case Swath::Method::RtkFloat:    return "RtkFloat";
    case Swath::Method::DR:          return "DR";
    case Swath::Method::Manual:      return "Manual";
    case Swath::Method::Sim:         return "Sim";
    case Swath::Method::PC:          return "PC";
    case Swath::Method::Other:       return "Other";
    default: return nullptr;
  }
} // Name(Swath::Method)

namespace tjg {
template<> struct EnumValues<Swath::Type>     : Swath::TypeList      { };
template<> struct EnumValues<Swath::Option>   : Swath::OptionList    { };
template<> struct EnumValues<Swath::Direction>: Swath::DirectionList { };
template<> struct EnumValues<Swath::Extension>: Swath::ExtensionList { };
template<> struct EnumValues<Swath::Method>   : Swath::MethodList    { };
} // tjg

Swath::Swath(const XmlNode& node)
  : id  {RequireAttr<std::string>(node, "A")}
  , type{RequireAttr<Type>(node, "C")}
{
  for (const auto& a: node.attributes()) {
    auto k = tjg::name(a);
    if (k == "A" || k == "C")
      continue;
    if      (k == "B") name      = tjg::get_attr<std::string>(a);
    else if (k == "D") option    = tjg::get_attr<Option>(a);
    else if (k == "E") direction = tjg::get_attr<Direction>(a);
    else if (k == "F") extension = tjg::get_attr<Extension>(a);
    else if (k == "G") heading   = tjg::get_attr<Angle>(a);
    else if (k == "H") radius    = tjg::get_attr<unsigned>(a);
    else if (k == "I") method    = tjg::get_attr<Method>(a);
    else if (k == "J") horizontalAccuracy = tjg::get_attr<double>(a);
    else if (k == "K") verticalAccuracy   = tjg::get_attr<double>(a);
    else if (k == "L") baseId    = tjg::get_attr<std::string>(a);
    else if (k == "M") srid      = tjg::get_attr<std::string>(a);
    else if (k == "N") leftRemaining  = tjg::get_attr<unsigned>(a);
    else if (k == "O") rightRemaining = tjg::get_attr<unsigned>(a);
    else otherAttrs.emplace_back(k, a.value());
  }
  for (const auto& c: node.children()) {
    if (c.type() != pugi::node_element)
      continue;
    auto k = tjg::name(c);
    if      (k == "LSG") paths.emplace_back(c);
    else if (k == "PLN") boundaries.emplace_back(c);
    else std::cerr << "Swath: ignored element: " << k << '\n';
  }
} // Swath::ctor

void Swath::dump(XmlNode& node) const {
  node.set_name("GPN");
  node.append_attribute("A") = id;
  if (!name.empty()) node.append_attribute("B") = name;
  node.append_attribute("C") = static_cast<int>(type);
  if (option)    node.append_attribute("D") = static_cast<int>(*option);
  if (direction) node.append_attribute("E") = static_cast<int>(*direction);
  if (extension) node.append_attribute("F") = static_cast<int>(*extension);
  if (heading)   node.append_attribute("G") = *heading;
  if (radius)    node.append_attribute("H") = *radius;
  if (method)    node.append_attribute("I") = static_cast<int>(*method);
  if (horizontalAccuracy)
      node.append_attribute("J") = *horizontalAccuracy;
  if (verticalAccuracy)
      node.append_attribute("K") = *verticalAccuracy;
  if (!baseId.empty()) node.append_attribute("L") = baseId;
  if (!srid.empty())   node.append_attribute("M") = srid;
  if (leftRemaining)
      node.append_attribute("N") = *leftRemaining;
  if (rightRemaining)
      node.append_attribute("O") = *rightRemaining;
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
  for (const auto& b: boundaries) {
    auto c = node.append_child("PLN");
    b.dump(c);
  }
  for (const auto& p: paths) {
    auto c = node.append_child("LSG");
    p.dump(c);
  }
} // Swath::dump

struct Guide {      // GGP
  std::string id;   // A
  std::string name; // B optional
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  std::vector<Swath> swaths;
  Guide() = default;
  explicit Guide(std::string_view id_, std::string_view name_="")
    : id{id_}, name{name_}
  {
    if (!std::regex_match(id, isoxml::GgpRe))
      throw std::invalid_argument{"Guide: invalid id: " + id};
  }
  Guide(const XmlNode& node);
  void dump(XmlNode& node) const;
}; // Guide

Guide::Guide(const XmlNode& node)
  : id  {RequireAttr<std::string>(node, "A")}
  , name{RequireAttr<std::string>(node, "B")}
{
  for (const auto& a: node.attributes()) {
    auto k = tjg::name(a);
    if (k == "A" || k == "B")
      continue;
    otherAttrs.emplace_back(k, a.value());
  }
  for (const auto& c: node.children()) {
    if (c.type() != pugi::node_element)
      continue;
    auto k = tjg::name(c);
    if (k == "GPN") {
      swaths.emplace_back(c);
      continue;
    }
    std::cerr << "Guide: ignored element: " << k << '\n';
  }
} // Guide::ctor

void Guide::dump(XmlNode& node) const {
  node.set_name("GGP");
  node.append_attribute("A") = id;
  if (!name.empty()) node.append_attribute("B") = name;
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
  for (const auto& s: swaths) {
    auto c = node.append_child("GPN");
    s.dump(c);
  }
} // Guide::dump

struct Field { // PFD
  std::string id;     // A
  std::string code;   // B optional
  std::string name;   // C
  unsigned area;      // D
  std::string ctrId;  // E optional
  std::string frmId;  // F optional
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  std::vector<Polygon> parts;
  std::vector<Guide> guides;
  Field() = default;
  Field(std::string_view id_, std::string_view name_="", unsigned area_=0)
    : id{id_}, name{name_}, area{area_}
  {
    if (!std::regex_match(id, isoxml::PfdRe))
      throw std::invalid_argument{"Field: invalid id: " + id};
  }
  explicit Field(const XmlNode& x);
  void dump(XmlNode& x) const;
}; // Field

Field::Field(const XmlNode& x)
  : id   {RequireAttr<std::string>(x, "A")}
  , name {RequireAttr<std::string>(x, "C")}
  , area {RequireAttr<unsigned>   (x, "D")}
{
  for (const auto& a: x.attributes()) {
    auto k = tjg::name(a);
    if (k == "A" || k == "C" || k == "D")
      continue;
    if      (k == "B") code  = tjg::get_attr<std::string>(a);
    else if (k == "E") ctrId = tjg::get_attr<std::string>(a);
    else if (k == "F") frmId = tjg::get_attr<std::string>(a);
    else otherAttrs.emplace_back(k, a.value());
  }
  for (const auto& c: x.children()) {
    if (c.type() != pugi::node_element)
      continue;
    auto k = tjg::name(c);
    if      (k == "PLN") parts.emplace_back(c);
    else if (k == "GGP") guides.emplace_back(c);
    else std::cerr << "Field: ignored element " << k << '\n';
  }
} // Field ctor

void Field::dump(XmlNode& node) const {
  node.set_name("PFD");
  node.append_attribute("A") = id;
  if (!code.empty()) node.append_attribute("B") = code;
  node.append_attribute("C") = name;
  node.append_attribute("D") = area;
  if (!ctrId.empty()) node.append_attribute("E") = ctrId;
  if (!frmId.empty()) node.append_attribute("F") = frmId;
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
  for (const auto& p: parts) {
    auto c = node.append_child("PLN");
    p.dump(c);
  }
  for (const auto& g: guides) {
    auto c = node.append_child("GGP");
    g.dump(c);
  }
} // Field::dump

struct Value { // VPN
  std::string id;
  int offset;
  double scale;
  int decimals;
  std::string units; // optional
  std::string color; // optional
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  Value() = default;
  explicit Value(const XmlNode& node);
  void dump(XmlNode& node) const;
}; // Value

Value::Value(const XmlNode& node)
  : id      {RequireAttr<std::string>(node, "A")}
  , offset  {RequireAttr<int        >(node, "B")}
  , scale   {RequireAttr<double     >(node, "C")}
  , decimals{RequireAttr<int        >(node, "D")}
{
  for (const auto& a: node.attributes()) {
    auto k = tjg::name(a);
    if (k == "A" || k == "B" || k == "C" || k == "D")
      continue;
    if      (k == "E") units = tjg::get_attr<std::string>(a);
    else if (k == "F") color = tjg::get_attr<std::string>(a);
    else otherAttrs.emplace_back(k, a.value());
  }
} // Value::ctor

void Value::dump(XmlNode& node) const {
  node.set_name("VPN");
  node.append_attribute("A") = id;
  node.append_attribute("B") = offset;
  node.append_attribute("C") = scale;
  node.append_attribute("D") = decimals;
  if (!units.empty()) node.append_attribute("E") = units;
  if (!color.empty()) node.append_attribute("F") = color;
} // Value::dump

struct FarmDb {
  int versionMajor       =  3;  // required
  int versionMinor       =  0;  // required
  int dataTransferOrigin = -1;  // optional, -1 means "unset"
  std::string swVendor;         // optional
  std::string swVersion;        // optional
  std::vector<std::pair<std::string,std::string>> otherAttrs;
  std::vector<Customer> customers;
  std::vector<Farm>     farms;
  std::vector<Field>    fields;
  std::vector<Value>    values;
  FarmDb() = default;
  explicit FarmDb(std::string_view version)
    : swVendor{"Ekman and Gober"}, swVersion{version}
    { }
  FarmDb(const XmlNode& node);
  void dump(XmlNode& node) const;
}; // FarmDb

FarmDb::FarmDb(const XmlNode& node)
  : versionMajor{RequireAttr<int>(node, isoxml::root_attr::VersionMajor)}
  , versionMinor{RequireAttr<int>(node, isoxml::root_attr::VersionMinor)}
{
  for (const auto& a : node.attributes()) {
    const auto k = tjg::name(a);
    if (   k == isoxml::root_attr::VersionMajor
        || k == isoxml::root_attr::VersionMinor)
      continue;
    if (k == isoxml::root_attr::DataTransferOrigin)
      dataTransferOrigin = tjg::get_attr<int>(a);
    else if (k == isoxml::root_attr::MgmtSoftwareManufacturer)
      swVendor   = tjg::get_attr<std::string>(a);
    else if (k == isoxml::root_attr::MgmtSoftwareVersion)
      swVersion        = tjg::get_attr<std::string>(a);
    else
      otherAttrs.emplace_back(k, a.value());
  }
  if (versionMajor < 0 || versionMinor < 0)
    throw std::runtime_error{"FarmDb: missing VersionMajor/VersionMinor"};
  for (const auto& c: node.children()) {
    if (c.type() != pugi::node_element)
      continue;
    auto k = tjg::name(c);
    if      (k == "CTR") customers.emplace_back(c);
    else if (k == "FRM") farms.emplace_back(c);
    else if (k == "PFD") fields.emplace_back(c);
    else if (k == "VPN") values.emplace_back(c);
    else std::cerr << "Root: ignored element: " << k << '\n';
  }
} // FarmDb ctor

void FarmDb::dump(XmlNode& node) const {
  node.set_name(isoxml::Root);
  if (versionMajor < 0 || versionMinor < 0) {
    auto msg = std::format("FarmDb::dump: invlid version: {}.{}",
                           versionMajor, versionMinor);
    throw std::runtime_error{msg};
  }
  node.append_attribute(isoxml::root_attr::VersionMajor) = versionMajor;
  node.append_attribute(isoxml::root_attr::VersionMinor) = versionMinor;
  if (dataTransferOrigin != -1) {
    node.append_attribute(isoxml::root_attr::DataTransferOrigin) =
                                                            dataTransferOrigin;
  }
  node.append_attribute(isoxml::root_attr::MgmtSoftwareManufacturer) =
                                                            swVendor;
  node.append_attribute(isoxml::root_attr::MgmtSoftwareVersion) =
                                                            swVersion;
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
  for (const auto& ctr: customers) {
    auto c = node.append_child("CTR");
    ctr.dump(c);
  }
  for (const auto& frm: farms) {
    auto c = node.append_child("FRM");
    frm.dump(c);
  }
  for (const auto& pfd: fields) {
    auto c = node.append_child("PFD");
    pfd.dump(c);
  }
  for (const auto& vpn: values) {
    auto c = node.append_child("VPN");
    vpn.dump(c);
  }
} // FarmDb::dump

// ---------------------------------------------------------------------
// Main: load, basic validation, ordered traversal with dispatch.
int main(int argc, const char *argv[]) {
  try {
    const auto path = (argc > 1) ? argv[1] : "TASKDATA.XML";

    auto doc = pugi::xml_document{};
    auto res = doc.load_file(path, pugi::parse_default | pugi::parse_ws_pcdata);
    if (!res) {
      std::cerr << std::format("XML parse error: {} (offset {})\n",
                               res.description(), res.offset);
      return 1;
    }

    auto root = doc.child(isoxml::Root);
    if (!root) {
      std::cerr << std::format("error: missing root <{}>\n", isoxml::Root);
      return 1;
    }

    auto top = FarmDb{root};
    std::cout << top.customers.size() << " customers\n"
              << top.farms.size()     << " farms\n"
              << top.fields.size()    << " fields\n"
              << top.values.size()    << " values\n";

    top.swVendor = "Terry Golubiewski";
    top.swVersion = "0.1 (alpha)";

    auto doc2 = pugi::xml_document{};
    auto root2 = doc2.append_child(isoxml::Root);
    top.dump(root2);

    auto ok = doc2.save_file("out.xml", "  ");
    if (!ok)
      throw std::runtime_error{"Error writing 'out.xml'"};
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

