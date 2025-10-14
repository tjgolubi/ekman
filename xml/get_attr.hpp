#include "enum_help.hpp"

#include <pugixml.hpp>

#include <string>
#include <string_view>
#include <optional>
#include <format>
#include <stdexcept>
#include <concepts>
#include <charconv>
#include <type_traits>
#include <utility>
#include <cctype>

namespace tjg {

inline std::string_view name(const pugi::xml_node& n) noexcept
  { return std::string_view{n.name()}; }
inline std::string_view value(const pugi::xml_node& n) noexcept
  { return std::string_view{n.value()}; }
inline std::string_view name(const pugi::xml_attribute& a) noexcept
  { return std::string_view{a.name()}; }
inline std::string_view value(const pugi::xml_attribute& a) noexcept
  { return std::string_view{a.value()}; }

namespace detail {

constexpr bool IsEqual(const char* a, const char* b) noexcept {
  if (!a || !b) return a == b;
  for (; *a && *b; ++a, ++b) {
    unsigned char ca = static_cast<unsigned char>(*a);
    unsigned char cb = static_cast<unsigned char>(*b);
    if (std::tolower(ca) != std::tolower(cb)) return false;
  }
  return (*a == *b);
} // IsEqual

template<class> inline constexpr bool DependentFalse = false;

template<class E>
concept HasFromChars = Enum<E> && requires(std::string_view sv) {
  { FromChars<E>(sv) } -> std::same_as<std::optional<E>>;
};

} // detail

template<typename T>
inline std::optional<T> try_get_attr(const pugi::xml_attribute& a) noexcept {
  if (!a) return std::nullopt;

  if constexpr (std::is_same_v<T, const char*>)
    { return a.value(); }
  else if constexpr (std::is_same_v<T, std::string_view>)
    { return std::string_view{a.value()}; }
  else if constexpr (std::is_same_v<T, std::string>)
    { return std::string{a.value()}; }
  else if constexpr (std::is_same_v<T, bool>) {
    auto s = a.value();
    if (!s) return std::nullopt;
    switch (s[0]) {
      case '1': if (!s[1]) return true;
                break;
      case '0': if (!s[1]) return false;
                break;
      case 't':
      case 'T': if (!s[1] || detail::IsEqual(s, "true")) return true;
                break;
      case 'f':
      case 'F': if (!s[1] || detail::IsEqual(s, "false")) return false;
                break;
      case 'y':
      case 'Y': if (!s[1] || detail::IsEqual(s, "yes")) return true;
                break;
      case 'n':
      case 'N': if (!s[1] || detail::IsEqual(s, "no")) return false;
                break;
      case 'o':
      case 'O': if (detail::IsEqual(s, "on" )) return true;
                if (detail::IsEqual(s, "off")) return false;
                break;
      default:  break;
    }
  }
  else if constexpr (Enum<T>) {
    if constexpr (detail::HasFromChars<T>) {
      auto v = FromChars<T>(std::string_view{a.value()});
      if (v) return *v;
    }
    if constexpr (EnumWithCast<T>) {
      using U = std::underlying_type_t<T>;
      if (auto u = try_get_attr<U>(a)) return enum_cast<T>(*u);
    }
  }
  else if constexpr (std::integral<T>) {
    auto s = a.value();
    if (!s || !s[0]) return std::nullopt;
    while (std::isspace(static_cast<unsigned char>(*s)))
      ++s;
    if (*s == '+') ++s;
    int base = 10;
    if constexpr (std::unsigned_integral<T>) {
      if (*s == '0') {
        switch (*++s) {
          case 'x': case 'X': ++s; base = 16; break;
          case 'b': case 'B': ++s; base =  2; break;
          case '\0': return T{0};
          default: base = 8; break;
        }
      }
    }
    auto sv = std::string_view{s};
    auto value = T{};
    const auto first = sv.data();
    const auto last  = first + sv.size();
    auto r = std::from_chars(first, last, value, base);
    if (r.ptr == last && r.ec == std::errc{})
      return value;
  }
  else if constexpr (std::floating_point<T>) {
    auto s = a.value();
    if (!s || !s[0]) return std::nullopt;
    while (std::isspace(static_cast<unsigned char>(*s)))
      ++s;
    if (*s == '+') ++s;
    auto sv = std::string_view{s};
    auto value = T{};
    const auto first = sv.data();
    const auto last  = first + sv.size();
    auto r = std::from_chars(first, last, value);
    if (r.ptr == last && r.ec == std::errc{})
      return value;
  }
  else {
    static_assert(detail::DependentFalse<T>, "get_attr: invalid type");
  }
  return std::nullopt;
} // try_get_attr

template<typename T>
inline T get_attr(const pugi::xml_attribute& a) {
  if (!a) throw std::runtime_error{"get_attr: empty attribute"};
  auto v = try_get_attr<T>(a);
  if (v)
    return *v;
  auto msg = std::format("get_attr: invalid attribute: {}={}",
                         a.name(), a.value());
  throw std::runtime_error{msg};
} // get_attr

} // tjg
