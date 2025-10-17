#pragma once

#include <numbers>
#include <compare>
#include <cmath>

namespace geom {

class Radians {
private:
  static constexpr auto _Pi = std::numbers::pi;
  static constexpr auto _TwoPi = 2.0 * _Pi;
  static constexpr auto _DegPerRad = 180.0 / _Pi;
  static constexpr auto _RadPerDeg = _Pi / 180.0;

  double _value = 0.0;

  constexpr Radians& _wrap(double x) noexcept {
    if (x > _Pi)
      x -= _TwoPi;
    else if (x <= -_Pi)
      x += _TwoPi;
    _value = x;
    return *this;
  }

  static constexpr double Wrap(double theta) noexcept {
    if (std::is_constant_evaluated()) {
      while (theta >   _Pi) theta -= _TwoPi;
      while (theta <= -_Pi) theta += _TwoPi;
      return theta;
    }
    else {
      return std::remainder(theta, _TwoPi);
    }
  } // Wrap

public:
  constexpr Radians() noexcept = default;
  constexpr explicit Radians(double theta) noexcept
    : _value{Wrap(theta)} { }
  static struct NoWrapT { } NoWrap;
  constexpr Radians(double v, NoWrapT) noexcept : _value{v} { }
  constexpr Radians(const Radians&) noexcept = default;
  constexpr Radians& operator=(const Radians&) noexcept = default;
  [[nodiscard]] constexpr bool operator==(const Radians&) const noexcept = default;
  [[nodiscard]] constexpr auto operator<=>(const Radians&) const noexcept = default;

  [[nodiscard]] static constexpr Radians FromDegrees(double deg) noexcept
    { return Radians{deg * _RadPerDeg}; }

  [[nodiscard]] constexpr double value() const noexcept { return _value; }
  [[nodiscard]] constexpr double degrees() const noexcept
    { return value() * _DegPerRad; }

  constexpr Radians& operator+=(const Radians& rhs) noexcept
    { return _wrap(_value + rhs._value); }
  constexpr Radians& operator-=(const Radians& rhs) noexcept
    { return _wrap(_value - rhs._value); }
  constexpr Radians& operator*=(double s) noexcept
    { *this = Radians{_value * s}; return *this; }
  constexpr Radians& operator/=(double s) noexcept
    { *this = Radians{_value / s}; return *this; }

  [[nodiscard]] constexpr Radians operator+() const noexcept { return *this; }
  [[nodiscard]] constexpr Radians operator-() const noexcept {
    if (_value == _Pi) [[unlikely]]
      return Radians{_Pi, NoWrap};
    return Radians{-_value, NoWrap};
  }

  [[nodiscard]] constexpr Radians operator+(const Radians& rhs) const noexcept
    { return Radians{*this} += rhs; }
  [[nodiscard]] constexpr Radians operator-(const Radians& rhs) const noexcept
    { return Radians{*this} -= rhs; }

  [[nodiscard]] constexpr Radians operator*(double s) const noexcept
    { return Radians{_value * s}; }
  [[nodiscard]] constexpr Radians operator/(double s) const noexcept
    { return Radians{_value / s}; }

}; // Radians

[[nodiscard]] constexpr Radians asin(double x) noexcept
  { return Radians{std::asin(x), Radians::NoWrap}; }
[[nodiscard]] constexpr Radians acos(double x) noexcept
  { return Radians{std::acos(x), Radians::NoWrap}; }
[[nodiscard]] constexpr Radians atan(double x) noexcept
  { return Radians{std::atan(x), Radians::NoWrap}; }
[[nodiscard]] constexpr Radians atan2(double y, double x) noexcept
  { return Radians{std::atan2(y, x), Radians::NoWrap}; }

[[nodiscard]] constexpr Radians operator*(double s, Radians rhs) noexcept
  { return rhs * s; }
[[nodiscard]] constexpr Radians abs(const Radians& x) noexcept
  { return (x.value() >= 0) ? x : -x; }

[[nodiscard]] constexpr double sin(const Radians& x) noexcept
  { return std::sin(x.value()); }
[[nodiscard]] constexpr double cos(const Radians& x) noexcept
  { return std::cos(x.value()); }
[[nodiscard]] constexpr double tan(const Radians& x) noexcept
  { return std::tan(x.value()); }

constexpr Radians Pi = Radians{std::numbers::pi};

} // geom
