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
inline bool is_group() {
    return Group<Set,Op>;
}

template <class Set, class Op>
concept AbelianGroup = Group<Set, Op> && Commutative<Op>;

template <class Set, class Op>
inline bool is_abelian_group() {
    return AbelianGroup<Set,Op>;
}

template <typename T>
class RealAddition {
public:
    using Domain = CartesianProduct<RealNumbers<T>,RealNumbers<T>>;
    using Target = RealNumbers<T>;

    static inline T identity() { return RealNumbers<T>::zero; }
    static consteval bool has_identity(const T& x) { 
        return apply(identity(), x) == x;
    }

    static inline T inverse(const T& x) { return -x; }
    static consteval bool has_inverse(const T& x) { 
        return apply(inverse(x),x) == identity();
    }

    static inline T apply(const T& x, const T& y) { return x+y; }
    static inline T apply(const Domain::type& x) { return x.first+x.second; }

    static consteval bool is_commutative(const T& x, const T& y) {
        return apply(x,y) == apply(y,x);
    }
    static constexpr bool is_associative(const T& x, const T& y, const T& z) {
        return apply(x,apply(y,z)) == apply(apply(x,y),z);
    }
};

template <typename T>
class RealMultiplication {
public:
    using Domain = CartesianProduct<RealNumbers<T>,RealNumbers<T>>;
    using Target = RealNumbers<T>;

    static inline T identity() { return RealNumbers<T>::one; }
    static consteval bool has_identity(const T& x) { 
        return apply(identity(), x) == x;
    }

    static inline T inverse (const T& x) {
        if (x == RealNumbers<T>::zero) {
            throw std::domain_error("Zero has no multiplicative inverse");
        }
        return identity/x;
    }
    static consteval bool has_inverse(const T& x) { 
        return apply(inverse(x),x) == identity();
    }

    static inline T apply(const T& x, const T& y) { return x*y; }
    static inline T apply(const Domain::type& x) { return x.first*x.second; }

    static consteval bool is_commutative(const T& x, const T& y) {
         return apply(x,y) == apply(y,x);
    }
    static constexpr bool is_associative(const T& x, const T& y, const T& z) { 
        return apply(x,apply(y,z)) == apply(apply(x,y),z); 
    }

    template <class Op>
    static constexpr bool distributes_through(const T& x, const T& y, const T& z) {
        return apply(z,Op::apply(x,y)) == Op::apply(apply(z,x), apply(z,y));
    }
};

