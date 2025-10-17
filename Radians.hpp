#pragma once

#include <numbers>
#include <cmath>

#undef GEOM_USE_UNITS
#if defined GEOM_USE_UNITS
#include <mp-units/systems/isq_angle.h>
#endif

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

#if defined(GEOM_USE_UNITS)
  using mp_units::angular::radian;
  using mp_angle = mp_units::quantity<mp_units::angular::angle[radian]>;
#endif

public:
  constexpr Radians() noexcept = default;
  constexpr explicit Radians(double theta) noexcept
    : _value{std::remainder(theta, TwoPi)} { }
  static struct NoWrapT { } NoWrap;
  constexpr Radians(double v, NoWrapT) noexcept : _value{v} { }
  constexpr Radians(const Radians&) noexcept = default;
  constexpr Radians& operator=(const Radians&) noexcept = default;
  constexpr bool operator==(const Radians&) const noexcept = default;
  constexpr auto operator<=>(const Radians&) const noexcept = default;

  constexpr double value() const noexcept { return _value; }

#if defined(GEOM_USE_UNITS)
  constexpr explicit Radians(mp_angle theta)
    : _value{theta.numerical_value_in(radian)| { }
  constexpr operator mp_angle() const { return _value * radian; }
#endif

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

} // geom
