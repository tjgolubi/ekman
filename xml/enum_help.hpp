#include <string_view>
#include <optional>
#include <concepts>
#include <type_traits>
#include <utility>

namespace tjg {

template<class E>
concept Enum = std::is_enum_v<E>;

template<class E> struct EnumValues; // specialize per enum

template<Enum E, E... Vs>
struct EnumList {
  using enum_type = E;
  static constexpr std::size_t size = sizeof...(Vs);
};

template<Enum E, E... Vs, class F>
constexpr void ForEachEnum(EnumList<E, Vs...>, F&& f) {
  (static_cast<void>(f(std::integral_constant<E, Vs>{})), ...);
}

template<class E>
concept EnumWithName = Enum<E> && requires(E e) {
  { Name(e) } -> std::same_as<const char*>;
};

template<EnumWithName E, typename T>
  requires std::convertible_to<T, std::underlying_type_t<E>>
constexpr std::optional<E> enum_cast(T&& other) noexcept {
  using U = std::underlying_type_t<E>;
  const auto u = static_cast<U>(std::forward<T>(other));
  const auto e = static_cast<E>(u);
  return Name(e) ? std::optional<E>{e} : std::nullopt;
} // enum_cast

template<class E>
concept EnumWithCast = Enum<E> && requires(std::underlying_type_t<E> x) {
  { enum_cast<E>(x) } -> std::same_as<std::optional<E>>;
};

namespace detail {

constexpr bool EqualSvCz(std::string_view sv, const char* cz) noexcept {
  if (!cz) return false;
  const char* p = cz;
  const auto n = sv.size();
  auto i = std::size_t{0};
  for ( ; i < n && *p; ++i, ++p) {
    if (sv[i] != *p) return false;
  }
  return (i == n) && (*p == '\0');
} // EqualSvCz

} // detail

template<Enum E>
  requires EnumWithName<E> && requires { typename EnumValues<E>; }
constexpr std::optional<E> FromChars(std::string_view s) {
  std::optional<E> out;
  ForEachEnum(EnumValues<E>{}, [&](auto V) {
    if (out) return;
    constexpr E e = V();
    if (const char* n = Name(e); n && detail::EqualSvCz(s, n))
      out = e;
  });
  return out;
} // FromChars

} // tjg
