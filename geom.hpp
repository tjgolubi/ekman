#pragma once
#include <boost/geometry.hpp>

#include <gsl-lite/gsl-lite.hpp>
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

  constexpr Radians& _wrap(double x) noexcept {
    if (x > Pi)
      x -= TwoPi;
    else if (x <= -Pi)
      x += TwoPi;
    _value = x;
    return *this;
  }

public:
  constexpr Radians() noexcept = default;
  constexpr explicit Radians(double x) noexcept : _value{std::remainder(x, TwoPi)} { }
  static struct NoWrapT { } NoWrap;
  constexpr Radians(double v, NoWrapT) noexcept : _value{v} { }
  constexpr Radians(const Radians&) noexcept = default;
  constexpr Radians& operator=(const Radians&) noexcept = default;
  constexpr bool operator==(const Radians&) const noexcept = default;
  constexpr auto operator<=>(const Radians&) const noexcept = default;

  constexpr double value() const noexcept { return _value; }

  constexpr Radians& operator+=(const Radians& rhs) noexcept
    { return _wrap(_value + rhs._value); }
  constexpr Radians& operator-=(const Radians& rhs) noexcept
    { return _wrap(_value - rhs._value); }
  constexpr Radians& operator*=(double s) noexcept
    { *this = Radians{_value * s}; return *this; }
  constexpr Radians& operator/=(double s) noexcept
    { *this = Radians{_value / s}; return *this; }

  constexpr Radians operator+() const noexcept { return *this; }
  constexpr Radians operator-() const noexcept {
    if (_value == Pi) [[unlikely]]
      return Radians{Pi, NoWrap};
    return Radians{-_value, NoWrap};
  }

  constexpr Radians operator+(const Radians& rhs) const noexcept
    { return Radians{*this} += rhs; }
  constexpr Radians operator-(const Radians& rhs) const noexcept
    { return Radians{*this} -= rhs; }

  constexpr Radians operator*(double s) const noexcept
    { return Radians{_value * s}; }
  constexpr Radians operator/(double s) const noexcept
    { return Radians{_value / s}; }

}; // Radians

constexpr Radians asin(double x) noexcept
  { return Radians{std::asin(x), Radians::NoWrap}; }
constexpr Radians acos(double x) noexcept
  { return Radians{std::acos(x), Radians::NoWrap}; }
constexpr Radians atan(double x) noexcept
  { return Radians{std::atan(x), Radians::NoWrap}; }
constexpr Radians atan2(double y, double x) noexcept
  { return Radians{std::atan2(y, x), Radians::NoWrap}; }

constexpr Radians operator*(double s, Radians rhs) noexcept { return rhs * s; }
constexpr Radians abs(const Radians& x) noexcept
  { return (x.value() >= 0) ? x : -x; }

constexpr double sin(const Radians& x) noexcept { return std::sin(x.value()); }
constexpr double cos(const Radians& x) noexcept { return std::cos(x.value()); }
constexpr double tan(const Radians& x) noexcept { return std::tan(x.value()); }

constexpr Radians FromDegrees(double deg) noexcept
    { return Radians{deg * RadPerDeg}; }

constexpr Radians ToRadians(double deg) noexcept { return FromDegrees(deg); }

constexpr double ToDegrees(Radians theta) noexcept
    { return theta.value() * DegPerRad; }

template<typename U> struct SquaredType { using type = U; };
template<typename U>
using SquaredTypeT = SquaredType<U>::type;

template<typename Units=double>
struct Vec {
  using value_type   = Units;
  using squared_type = SquaredTypeT<Units>;
  value_type dx = value_type{};
  value_type dy = value_type{};
  constexpr Vec() noexcept = default;
  constexpr Vec(const Vec&) noexcept = default;
  constexpr Vec& operator=(const Vec&) noexcept = default;
  constexpr bool operator==(const Vec&) const noexcept = default;
  constexpr Vec(value_type dx_, value_type dy_) : dx{dx_}, dy{dy_} { }
  constexpr Vec(value_type mag, Radians theta)
    : dx{mag * cos(theta)}, dy{mag * sin(theta)}
    { gsl_Expects(mag >= value_type{0}); }
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
  constexpr squared_type norm2() const noexcept { return dx*dx + dy*dy; }
  constexpr value_type norm() const noexcept { return std::hypot(dx, dy); }
  constexpr Vec unit() const noexcept { return *this / norm(); }
  constexpr Radians angle() const noexcept
    { return geom::atan2(dy, dx); }
  friend constexpr Vec operator*(double s, const Vec& rhs)
    { return rhs * s; }
  friend constexpr squared_type dot(const Vec& u, const Vec& v)
    { return u.dx * v.dx + u.dy * v.dy; }
  friend constexpr squared_type operator*(const Vec& u, const Vec& v)
    { return dot(u, v); }
  friend constexpr squared_type cross(const Vec& u, const Vec& v) noexcept
    { return u.dx * v.dy - u.dy * v.dx; }
  constexpr Radians angle_wrt(const Vec& ref) const noexcept
    { return geom::atan2(cross(ref, *this), dot(ref, *this)); }
}; // Vec

template<typename Units=double>
struct Pt {
  Units x = Units{};
  Units y = Units{};
  constexpr Pt() = default;
  constexpr Pt(const Pt&) = default;
  constexpr Pt& operator=(const Pt&) = default;
  constexpr bool operator==(const Pt&) const = default;
  constexpr Pt(Units x_, Units y_) : x{x_}, y{y_} { }
}; // Pt

} // geom

namespace boost::geometry::traits {

template<typename U>
struct tag<geom::Pt<U>> { using type = point_tag; };

template<typename U>
struct coordinate_type<geom::Pt<U>> { using type = U; };

template<typename U>
struct coordinate_system<geom::Pt<U>> { using type = cs::cartesian; };

template<typename U>
struct dimension<geom::Pt<U>> : std::integral_constant<std::size_t, 2> { };

template<std::size_t Dim, typename U>
requires (Dim == 0 || Dim == 1)
struct access<geom::Pt<U>, Dim> {
  static constexpr U get(const geom::Pt<U>& p) {
    if constexpr (Dim == 0)
      return p.x;
    else
      return p.y;
  }
  static constexpr void set(geom::Pt<U>& p, U v) {
    if constexpr (Dim == 0)
      p.x = v;
    else
      p.y = v;
  }
}; // access

template<typename U>
struct make<geom::Pt<U>> {
  using point_type = geom::Pt<U>;
  static constexpr auto is_specialized = true;
  static constexpr point_type apply(U x, U y)
    { return point_type{x, y}; }
}; // make

} // boost::geometry::traits

namespace geom {

template<typename U>
constexpr Vec<U> operator-(const Pt<U>& lhs, const Pt<U>& rhs) noexcept
  { return Vec<U>{lhs.x-rhs.x, lhs.y-rhs.y}; }

template<typename U>
constexpr Pt<U> operator+(const Pt<U>& p, const Vec<U>& v) noexcept
  { return Pt<U>{p.x + v.dx, p.y + v.dy}; }

template<typename U>
constexpr Pt<U> operator-(const Pt<U>& p, const Vec<U>& v) noexcept
  { return Pt<U>{p.x - v.dx, p.y - v.dy}; }

template<typename U>
constexpr U Dist(const Pt<U>& a, const Pt<U>& b) noexcept
  { return (b - a).norm(); }

template<typename U>
constexpr SquaredTypeT<U> Dist2(const Pt<U>& a, const Pt<U>& b) noexcept
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
static_assert(vh.angle() == geom::atan2(4.0, 3.0));
static_assert(vh.angle_wrt(vx) == vh.angle());
static_assert(vh.angle_wrt(vy) == vh.angle() - Radians{HalfPi});
#endif
} // test

} // geom
