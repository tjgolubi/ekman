#pragma once
#include <boost/geometry.hpp>

#include <vector>
#include <numbers>
#include <cmath>

namespace geom {

constexpr auto Pi = std::numbers::pi;
constexpr auto TwoPi = 2.0 * Pi;
constexpr auto HalfPi = Pi / 2.0;
constexpr auto DegPerRad = 180.0 / Pi;
constexpr auto RadPerDeg = Pi / 180.0;

class Radians {
  double _value = 0.0;

public:
  constexpr Radians() noexcept = default;
  constexpr explicit Radians(double v) noexcept : _value{v} { }
  constexpr Radians(const Radians&) noexcept = default;
  constexpr Radians& operator=(const Radians&) noexcept = default;
  constexpr bool operator==(const Radians&) const noexcept = default;

  friend constexpr Radians Normalize(Radians x) noexcept {
    if (std::is_constant_evaluated()) {
      while (x._value >   Pi) x._value -= TwoPi;
      while (x._value <= -Pi) x._value += TwoPi;
    }
    else {
      x._value = std::fmod(x._value + Pi, TwoPi);
      x._value += (x._value < 0.0) ? Pi : -Pi;
    }
    return x;
  } // Normalize

  constexpr Radians& normalize() noexcept { return *this = Normalize(*this); }
  constexpr Radians  normalize() const noexcept { return Normalize(*this); }

  constexpr double value() const noexcept { return normalize()._value; }
  constexpr double value()       noexcept { return normalize()._value; }

  constexpr Radians& operator+=(const Radians& rhs) noexcept
    { _value += rhs._value; return *this; }
  constexpr Radians& operator-=(const Radians& rhs) noexcept
    { _value -= rhs._value; return *this; }

  constexpr Radians operator+() const noexcept { return *this; }
  constexpr Radians operator-() const noexcept { return Radians{-_value}; }

  constexpr Radians operator+(const Radians& rhs) const noexcept
    { return Radians{_value + rhs._value}; }
  constexpr Radians operator-(const Radians& rhs) const noexcept
    { return Radians{_value - rhs._value}; }

  constexpr Radians operator*(double s) const noexcept
    { return Radians{_value * s}; }
  constexpr Radians operator/(double s) const noexcept
    { return Radians{_value / s}; }
  constexpr Radians& operator*=(double s) noexcept
    { _value *= s; return *this; }
  constexpr Radians& operator/=(double s) noexcept
    { _value /= s; return *this; }

  friend constexpr Radians operator*(double s, Radians rhs) noexcept
    { return rhs * s; }
  friend constexpr double sin(Radians x) noexcept { return std::sin(x._value); }
  friend constexpr double cos(Radians x) noexcept { return std::cos(x._value); }
  friend constexpr double tan(Radians x) noexcept { return std::tan(x._value); }
  friend constexpr Radians abs(Radians x) noexcept
    { x.normalize(); return (x._value >= 0) ? x : -x; }
}; // Radians

constexpr Radians FromDegrees(double deg) noexcept
    { return Radians{deg * RadPerDeg}; }

constexpr Radians ToRadians(double deg) noexcept { return FromDegrees(deg); }

constexpr double ToDegrees(Radians theta) noexcept
    { return theta.value() * DegPerRad; }

struct Vec {
  double dx=0.0;
  double dy=0.0;
  constexpr Vec() noexcept = default;
  constexpr Vec(const Vec&) noexcept = default;
  constexpr Vec& operator=(const Vec&) noexcept = default;
  constexpr bool operator==(const Vec&) const noexcept = default;
  constexpr Vec(double dx_, double dy_) : dx{dx_}, dy{dy_} { }
  constexpr Vec(double mag, Radians theta)
    : dx{mag * cos(theta)}, dy{mag * sin(theta)} { }
  constexpr const Vec& operator+() const noexcept { return *this; }
  constexpr Vec operator-() const noexcept { return Vec{-dx, -dy}; }
  constexpr Vec& operator+=(const Vec& rhs) noexcept
    { dx += rhs.dx; dy += rhs.dy; return *this; }
  constexpr Vec& operator-=(const Vec& rhs) noexcept
    { dx -= rhs.dx; dy -= rhs.dy; return *this; }
  constexpr Vec& operator*=(double s) noexcept
    { dx *= s; dy *= s; return *this; }
  constexpr Vec& operator/=(double s) noexcept
    { dx /= s; dy /= s; return *this; }
  constexpr Vec operator+(const Vec& rhs) const noexcept
    { return Vec{dx+rhs.dx, dy+rhs.dy}; }
  constexpr Vec operator-(const Vec& rhs) const noexcept
    { return Vec{dx-rhs.dx, dy-rhs.dy}; }
  constexpr Vec operator*(double s) const noexcept
    { return Vec{dx*s, dy*s}; }
  constexpr Vec operator/(double s) const noexcept
    { return Vec{dx/s, dy/s}; }
  constexpr double norm2() const noexcept { return dx*dx + dy*dy; }
  constexpr double norm() const noexcept { return std::hypot(dx, dy); }
  constexpr Vec unit() const noexcept { return *this / norm(); }
  constexpr Radians angle() const noexcept
    { return Radians{std::atan2(dy, dx)}; }
  friend constexpr Vec operator*(double s, const Vec& rhs)
    { return rhs * s; }
  friend constexpr double dot(const Vec& u, const Vec& v)
    { return u.dx * v.dx + u.dy * v.dy; }
  friend constexpr double operator*(const Vec& u, const Vec& v)
    { return dot(u, v); }
  friend constexpr double cross(const Vec& u, const Vec& v) noexcept
    { return u.dx * v.dy - u.dy * v.dx; }
  constexpr Radians angle_wrt(const Vec& ref) const noexcept
    { return Radians{std::atan2(cross(ref, *this), dot(ref, *this))}; }
}; // Vec

struct Pt {
  double x = 0.0;
  double y = 0.0;
  constexpr Pt() = default;
  constexpr Pt(const Pt&) = default;
  constexpr Pt& operator=(const Pt&) = default;
  constexpr bool operator==(const Pt&) const = default;
  constexpr Pt(double x_, double y_) : x{x_}, y{y_} { }
}; // Pt

} // geom

namespace boost::geometry::traits {

template<>
struct tag<geom::Pt> { using type = point_tag; };

template<>
struct coordinate_type<geom::Pt> { using type = double; };

template<>
struct coordinate_system<geom::Pt> { using type = cs::cartesian; };

template<>
struct dimension<geom::Pt> : std::integral_constant<std::size_t, 2> { };

template<std::size_t Dim>
requires (Dim == 0 || Dim == 1)
struct access<geom::Pt, Dim> {
  static constexpr double get(const geom::Pt& p) {
    if constexpr (Dim == 0)
      return p.x;
    else
      return p.y;
  }
  static constexpr void set(geom::Pt& p, double v) {
    if constexpr (Dim == 0)
      p.x = v;
    else
      p.y = v;
  }
}; // access

template<>
struct make<geom::Pt> {
  using point_type = geom::Pt;
  static constexpr auto is_specialized = true;
  static constexpr point_type apply(double x, double y)
    { return point_type{x, y}; }
}; // make

} // boost::geometry::traits

namespace geom {

constexpr Vec operator-(const Pt& lhs, const Pt& rhs) noexcept
  { return Vec{lhs.x-rhs.x, lhs.y-rhs.y}; }

constexpr Pt operator+(const Pt& p, const Vec& v) noexcept
  { return Pt{p.x + v.dx, p.y + v.dy}; }

constexpr Pt operator-(const Pt& p, const Vec& v) noexcept
  { return Pt{p.x - v.dx, p.y - v.dy}; }

constexpr double Dist(const Pt& a, const Pt& b) noexcept
  { return (b - a).norm(); }

constexpr double Dist2(const Pt& a, const Pt& b) noexcept
  { return (b - a).norm2(); }

namespace test {
constexpr auto p0 = Pt{0, 0};
constexpr auto p1 = Pt{3, 0};
constexpr auto p2 = Pt{3, 4};
constexpr auto vx = p1 - p0;
constexpr auto vy = p2 - p1;
constexpr auto vh = vx + vy;
static_assert(p0 + vh == p2);
static_assert(vh.norm2() == 25.0);
#if 0
static_assert(vh.norm() == 5.0);
static_assert(vx.angle() == Radians{0.0});
static_assert(vy.angle() == Radians{HalfPi});
static_assert(vh.angle() == Radians{std::atan2(4.0, 3.0)});
static_assert(vh.angle_wrt(vx) == vh.angle());
static_assert(vh.angle_wrt(vy) == vh.angle() - Radians{HalfPi});
#endif
} // test

} // geom
