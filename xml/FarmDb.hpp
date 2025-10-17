/// @file
/// @copyright 2025 Terry Golubiewski, all rights reserved.
/// @author Terry Golubiewski

#pragma once

#include "enum_help.hpp"

#include <string>
#include <vector>
#include <string_view>
#include <span>
#include <optional>
#include <functional>
#include <utility>

namespace pugi {
struct xml_node;
struct xml_attribute;
} // pugi

namespace farm_db {

using XmlNode = pugi::xml_node;
using XmlAttr = pugi::xml_attribute;

struct Customer {    // CTR
  std::string id;    // A
  std::string name;  // B
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  Customer() = default;
  Customer(std::string_view id_, std::string_view name_);
  explicit Customer(const XmlNode& node);
  void dump(XmlNode& node) const;
}; // Customer

struct Farm {           // FRM
  std::string id;       // A
  std::string name;     // B
  std::string ctrId;    // I optional
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  Farm() = default;
  Farm(std::string_view id_, std::string_view name_);
  explicit Farm(const XmlNode& node);
  void dump(XmlNode& node) const;
}; // Farm

using Degrees = double;
using Angle   = Degrees;

struct LatLon {
  Angle latitude;
  Angle longitude;
  LatLon(Angle lat_, Angle lon_) : latitude{lat_}, longitude(lon_) { }
}; // LatLon

using Path = std::vector<LatLon>;

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
  Angle latitude()  const { return point.latitude;  }
  Angle longitude() const { return point.longitude; }
  Point() = default;
  Point(const LatLon& pt, Type type_ = Type::Other)
    : type{type_}, point{pt} { }
  Point(Angle lat, Angle lon, Type type_ = Type::Other)
    : type{type_}, point{lon, lat} { }
  explicit Point(const XmlNode& x, TypeValidator validator=AnyType);
  void dump(XmlNode& x) const;
}; // Point

const char* Name(Point::Type x) noexcept;

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
  LineString(Type type_, Point::Type ptType, const Path& pts);

  explicit LineString(const XmlNode& x,
                      Point::TypeValidator point_validator=Point::AnyType);
  void dump(XmlNode& node) const;
}; // LineString

const char* Name(LineString::Type x) noexcept;

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
  explicit Polygon(const Path& poly, Type type_ = Type::Boundary);
  Polygon(const Path& poly, std::span<Path> inners_, Type type_ = Type::Boundary);
  explicit Polygon(const XmlNode& x);
  void dump(XmlNode& node) const;
}; // Polygon

const char* Name(Polygon::Type x) noexcept;

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
  void push(const Path& path);
  Swath() = default;
  explicit Swath(std::string_view id_, Type type_ = Type::Curve);
  Swath(std::string_view id_, const Path& ls, Type type_=Type::Curve);
  Swath(std::string_view id_, const std::span<Path> lines, Type type_=Type::Curve);
  explicit Swath(const XmlNode& node);
  void dump(XmlNode& node) const;
}; // Swath

const char* Name(Swath::Type x) noexcept;
const char* Name(Swath::Option x) noexcept;
const char* Name(Swath::Direction x) noexcept;
const char* Name(Swath::Extension x) noexcept;
const char* Name(Swath::Method x) noexcept;

struct Guide {      // GGP
  std::string id;   // A
  std::string name; // B optional
  std::vector<std::pair<std::string, std::string>> otherAttrs;
  std::vector<Swath> swaths;
  Guide() = default;
  explicit Guide(std::string_view id_, std::string_view name_="");
  Guide(const XmlNode& node);
  void dump(XmlNode& node) const;
}; // Guide

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
  Field(std::string_view id_, std::string_view name_="", unsigned area_=0);
  explicit Field(const XmlNode& x);
  void dump(XmlNode& x) const;
}; // Field

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
  explicit FarmDb(std::string_view version);
  FarmDb(const XmlNode& node);
  void dump(XmlNode& node) const;
}; // FarmDb

} // farm_db

namespace tjg {
template<>
struct EnumValues<farm_db::Point::Type>     : farm_db::Point::Types         { };
template<>
struct EnumValues<farm_db::LineString::Type>: farm_db::LineString::Types    { };
template<>
struct EnumValues<farm_db::Polygon::Type>   : farm_db::Polygon::Types       { };
template<>
struct EnumValues<farm_db::Swath::Type>     : farm_db::Swath::TypeList      { };
template<>
struct EnumValues<farm_db::Swath::Option>   : farm_db::Swath::OptionList    { };
template<>
struct EnumValues<farm_db::Swath::Direction>: farm_db::Swath::DirectionList { };
template<>
struct EnumValues<farm_db::Swath::Extension>: farm_db::Swath::ExtensionList { };
template<>
struct EnumValues<farm_db::Swath::Method>   : farm_db::Swath::MethodList    { };
} // tjg
