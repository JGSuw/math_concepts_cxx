#pragma once
#include <exception>
#include "operators.hpp"

template <class Set, class Op>
concept Group = PointSet<Set> && 
        Associative<Op> &&
        HasIdentityElement<Op> &&
        HasInverseElement<Op> &&
        MapsTo<Op,typename Op::Domain, Set>;

template <class Set, class Op>
consteval bool is_group() {
    return Group<Set,Op>;
}

template <class Set, class Op>
concept AbelianGroup = Group<Set, Op> && Commutative<Op>;

template <class Set, class Op>
consteval bool is_abelian_group() {
    return AbelianGroup<Set,Op>;
}

template <typename T>
class RealAddition {
public:
    using Domain = CartesianProduct<RealNumbers<T>,RealNumbers<T>>;
    using Target = RealNumbers<T>;
    static consteval T identity() { return RealNumbers<T>::zero; }
    static consteval T inverse(const T& x) { return -x; }
    static consteval T apply(const Domain::type& x) { return x.first+x.second; }
    static constexpr bool is_commutative = true;
    static constexpr bool is_associative = true;
};

template <typename T>
class RealMultiplication {
public:
    using Domain = CartesianProduct<RealNumbers<T>,RealNumbers<T>>;
    using Target = RealNumbers<T>;
    static consteval T identity() { return RealNumbers<T>::one; }
    static consteval T inverse (const T& x) {
        if (x == RealNumbers<T>::zero) {
            throw std::domain_error("Zero has no multiplicative inverse");
        }
        return identity/x;
    }
    static consteval T apply(const Domain::type& x) { return x.first*x.second; }
    static constexpr bool is_commutative = true;
    static constexpr bool is_associative = true;
    static constexpr bool distributes_through(RealAddition<T> addition) { return true; }
};

