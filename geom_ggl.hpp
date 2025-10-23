/// @file
/// @copyright 2025 Terry Golubiewski, all rights reserved.
/// @author Terry Golubiewski

#pragma once
#include "geom.hpp"

#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/core/cs.hpp>

#include <type_traits>
#include <cstddef>

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
  static constexpr U get(const geom::Pt<U>& p) { return get<Dim>(p); }
  static constexpr void set(geom::Pt<U>& p, U v) { get<Dim>(p) = v; }
}; // access

template<typename U>
struct make<geom::Pt<U>> {
  using point_type = geom::Pt<U>;
  static constexpr auto is_specialized = true;
  static constexpr point_type apply(U x, U y) { return point_type{x, y}; }
}; // make

} // boost::geometry::traits
