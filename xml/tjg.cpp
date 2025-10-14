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

#include <unordered_map>
#include <string>
#include <vector>
#include <string_view>
#include <functional>
#include <format>
#include <iostream>
#include <iomanip>
#include <utility>
#include <cstdint>

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
  explicit Customer(const XmlNode& node);
  void dump(XmlNode node) const;
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

void Customer::dump(XmlNode node) const {
  node.set_name("CTR");
  node.append_attribute("A") = id;
  node.append_attribute("B") = name;
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
} // Customer::dump

struct Farm {
  std::string id;
  std::string name;
  std::string ctrId;
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  Farm() = default;
  explicit Farm(const XmlNode& node);
  void dump(XmlNode node) const;
}; // Farm

Farm::Farm(const XmlNode& node)
  : id   {RequireAttr<std::string>(node, "A")}
  , name {RequireAttr<std::string>(node, "B")}
  , ctrId{RequireAttr<std::string>(node, "I")}
{
  for (const auto& a: node.attributes()) {
    auto k = tjg::name(a);
    if (k == "A" || k == "B" || k == "C")
      continue;
    otherAttrs.emplace_back(a.name(), a.value());
  }
} // Farm::ctor

void Farm::dump(XmlNode node) const {
  node.set_name("FRM");
  node.append_attribute("A") = id;
  node.append_attribute("B") = name;
  node.append_attribute("I") = ctrId;
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
} // Farm::dump

using Angle = double;

struct LatLon {
  Angle latitude  = Angle{};
  Angle longitude = Angle{};
  LatLon() = default;
  LatLon(Angle lat, Angle lon);
}; // LatLon

LatLon::LatLon(Angle lat, Angle lon)
  : latitude{lat}, longitude{lon}
{
  if ( lat <  -90.0 || lat >  +90.0
    || lon < -180.0 || lon > +180.0)
  {
    auto msg = std::format("LatLong: invalid ({}, {})", lat, lon);
    throw std::range_error{msg};
  }
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
  LatLon point;
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  using TypeValidator = std::function<bool(Type)>;
  static constexpr bool AnyType(Type) { return true; }
  Point() = default;
  explicit Point(const XmlNode& x, TypeValidator validator=AnyType);
  void dump(XmlNode x) const;
}; // Point

constexpr const char* Name(Point::Type x) noexcept {
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
  , point{RequireAttr<Angle>(x, "C"),
          RequireAttr<Angle>(x, "D")}
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

void Point::dump(XmlNode x) const {
  x.set_name("PNT");
  x.append_attribute("A") = static_cast<int>(type);
  x.append_attribute("C") = point.latitude;
  x.append_attribute("D") = point.longitude;
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

  LineString() = default;
  explicit LineString(const XmlNode& x,
                      Point::TypeValidator point_validator=Point::AnyType);
  void dump(XmlNode node) const;
}; // LineString

constexpr const char* Name(LineString::Type x) noexcept {
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

void LineString::dump(XmlNode node) const {
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

  Polygon() = default;
  explicit Polygon(const XmlNode& x);
  void dump(XmlNode node) const;
}; // Polygon

constexpr const char* Name(Polygon::Type x) noexcept {
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

void Polygon::dump(XmlNode node) const {
  node.set_name("PLN");
  node.append_attribute("A") = static_cast<int>(type);
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
  outer.dump(node.append_child("LSG"));
  for (const auto& lsg: inners)
    lsg.dump(node.append_child("LSG"));
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
  std::string name;
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
  Swath() = default;
  explicit Swath(const XmlNode& node);
  void dump(XmlNode node) const;
}; // Swath

constexpr const char* Name(Swath::Type x) noexcept {
  switch (x) {
    case Swath::Type::AB:     return "AB";
    case Swath::Type::APlus:  return "APlus";
    case Swath::Type::Curve:  return "Curve";
    case Swath::Type::Pivot:  return "Pivot";
    case Swath::Type::Spiral: return "Spiral";
    default: return nullptr;
  }
} // Name(Swath::Type)

constexpr const char* Name(Swath::Option x) noexcept {
  switch (x) {
    case Swath::Option::CW:   return "CW";
    case Swath::Option::CCW:  return "CCW";
    case Swath::Option::Full: return "Full";
    default: return nullptr;
  }
} // Name(Swath::Option)

constexpr const char* Name(Swath::Direction x) noexcept {
  switch (x) {
    case Swath::Direction::Both:  return "Both";
    case Swath::Direction::Left:  return "Left";
    case Swath::Direction::Right: return "Right";
    case Swath::Direction::None:  return "None";
    default: return nullptr;
  }
} // Name(Swath::Direction)

constexpr const char* Name(Swath::Extension x) noexcept {
  switch (x) {
    case Swath::Extension::Both:  return "Both";
    case Swath::Extension::First: return "First";
    case Swath::Extension::Last:  return "Last";
    case Swath::Extension::None:  return "None";
    default: return nullptr;
  }
} // Name(Swath::Extension)

constexpr const char* Name(Swath::Method x) noexcept {
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
  , name{RequireAttr<std::string>(node, "B")}
  , type{RequireAttr<Type>(node, "C")}
{
  for (const auto& a: node.attributes()) {
    auto k = tjg::name(a);
    if (k == "A" || k == "B" || k == "C")
      continue;
    if      (k == "D") option    = tjg::get_attr<Option>(a);
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

void Swath::dump(XmlNode node) const {
  node.set_name("GPN");
  node.append_attribute("A") = id;
  node.append_attribute("B") = name;
  node.append_attribute("C") = static_cast<int>(type);
  if (option)    node.append_attribute("D") = static_cast<int>(*option);
  if (direction) node.append_attribute("E") = static_cast<int>(*direction);
  if (extension) node.append_attribute("F") = static_cast<int>(*extension);
  if (heading)   node.append_attribute("G") = static_cast<double>(*heading);
  if (radius)    node.append_attribute("H") = static_cast<unsigned>(*radius);
  if (method)    node.append_attribute("I") = static_cast<int>(*method);
  if (horizontalAccuracy)
      node.append_attribute("J") = static_cast<double>(*horizontalAccuracy);
  if (verticalAccuracy)
      node.append_attribute("K") = static_cast<double>(*verticalAccuracy);
  if (!baseId.empty()) node.append_attribute("L") = baseId;
  if (!srid.empty())   node.append_attribute("M") = srid;
  if (leftRemaining)
      node.append_attribute("N") = static_cast<unsigned>(*leftRemaining);
  if (rightRemaining)
      node.append_attribute("O") = static_cast<unsigned>(*rightRemaining);
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
  for (const auto& b: boundaries)
    b.dump(node.append_child("PLN"));
  for (const auto& p: paths)
    p.dump(node.append_child("LSG"));
} // Swath::dump

struct Guide {      // GGP
  std::string id;   // A
  std::string name; // B
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  std::vector<Swath> swaths;
  Guide() = default;
  Guide(const XmlNode& node);
  void dump(XmlNode node) const;
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

void Guide::dump(XmlNode node) const {
  node.set_name("GGP");
  node.append_attribute("A") = id;
  node.append_attribute("B") = name;
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
  for (const auto& s: swaths)
    s.dump(node.append_child("GPN"));
} // Guide::dump

struct Field { // PFD
  std::string id;     // A
  std::string name;   // C
  unsigned area;      // D
  std::string ctrId;  // E
  std::string frmId;  // F
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  std::vector<Polygon> parts;
  std::vector<Guide> guides;
  Field() = default;
  explicit Field(const XmlNode& x);
  void dump(XmlNode x) const;
}; // Field

Field::Field(const XmlNode& x)
  : id   {RequireAttr<std::string>(x, "A")}
  , name {RequireAttr<std::string>(x, "C")}
  , area {RequireAttr<unsigned>   (x, "D")}
  , ctrId{RequireAttr<std::string>(x, "E")}
  , frmId{RequireAttr<std::string>(x, "F")}
{
  for (const auto& a: x.attributes()) {
    auto k = tjg::name(a);
    if (k == "A" || k == "C" || k == "D" || k == "E" || k == "F")
      continue;
    otherAttrs.emplace_back(k, a.value());
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

void Field::dump(XmlNode node) const {
  node.set_name("PFD");
  node.append_attribute("A") = id;
  node.append_attribute("C") = name;
  node.append_attribute("D") = area;
  node.append_attribute("E") = ctrId;
  node.append_attribute("F") = frmId;
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
  for (const auto& p: parts)
    p.dump(node.append_child("PLN"));
  for (const auto& g: guides)
    g.dump(node.append_child("GGP"));
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
  void dump(XmlNode node) const;
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

void Value::dump(XmlNode node) const {
  node.set_name("VPN");
  node.append_attribute("A") = id;
  node.append_attribute("B") = offset;
  node.append_attribute("C") = scale;
  node.append_attribute("D") = decimals;
  if (!units.empty()) node.append_attribute("E") = units;
  if (!color.empty()) node.append_attribute("F") = color;
} // Value::dump

struct RootMeta {
  int versionMajor       = -1;  // required
  int versionMinor       = -1;  // required
  int dataTransferOrigin = -1;  // optional, -1 means "unset"
  std::string mgmtManufacturer; // optional
  std::string mgmtVersion;      // optional
  std::vector<std::pair<std::string,std::string>> otherAttrs;
  std::vector<Customer> customers;
  std::vector<Farm>     farms;
  std::vector<Field>    fields;
  std::vector<Value>    values;
  RootMeta() = default;
  RootMeta(const XmlNode& node);
  void dump(XmlNode node) const;
}; // RootMeta

RootMeta::RootMeta(const XmlNode& node)
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
      mgmtManufacturer   = tjg::get_attr<std::string>(a);
    else if (k == isoxml::root_attr::MgmtSoftwareVersion)
      mgmtVersion        = tjg::get_attr<std::string>(a);
    else
      otherAttrs.emplace_back(k, a.value());
  }
  if (versionMajor < 0 || versionMinor < 0)
    throw std::runtime_error{"RootMeta: missing VersionMajor/VersionMinor"};
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
} // RootMeta ctor

void RootMeta::dump(XmlNode node) const {
  if (versionMajor < 0 || versionMinor < 0) {
    auto msg = std::format("RootMeta::dump: invlid version: {}.{}",
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
                                                            mgmtManufacturer;
  node.append_attribute(isoxml::root_attr::MgmtSoftwareVersion) =
                                                            mgmtVersion;
  for (const auto& [k, v]: otherAttrs)
    node.append_attribute(k) = v;
  for (const auto& ctr: customers)
    ctr.dump(node.append_child("CTR"));
  for (const auto& frm: farms)
    frm.dump(node.append_child("FRM"));
  for (const auto& pfd: fields)
    pfd.dump(node.append_child("PFD"));
} // RootMeta::dump

// ---------------------------------------------------------------------
// Main: load, basic validation, ordered traversal with dispatch.
int main(int argc, const char *argv[]) {
  const auto path = (argc > 1) ? argv[1] : "TASKDATA.XML";

  pugi::xml_document doc;
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

  auto input = RootMeta{root};
  std::cout << input.customers.size() << " customers\n"
            << input.farms.size()     << " farms\n"
            << input.fields.size()    << " fields\n"
            << input.values.size()    << " values\n";

  return 0;
} // main

