#pragma once
#include "fields.hpp"
#include <Eigen/Dense>


template <class A>
concept VectorSpace = 
        Field<typename A::F> &&
        AbelianGroup<typename A::V, typename A::Add> &&
        MapsTo<
            typename A::Mul,
            CartesianProduct<typename A::F::Set, typename A::V>, 
            typename A::V> &&
requires (A::F::Set::type& a, A::F::Set::type& b, A::V::type& x, A::V::type& y) {
    A::Mul::distributes_through(a,x,y);
    A::Mul::distributes_through(a,b,x);
};

template <class A>
bool is_vectorspace() {
    return VectorSpace<A>;
};

template <class A, class B>
concept Dual = 
        VectorSpace<A> &&
        VectorSpace<B> &&
        std::same_as<typename A::F, typename B::F> &&
requires(A::F::Set::type& x, A::V::type& covector, B::V::type& vector) {
    x = A::Contraction::apply(covector, vector);
};

template <VectorSpace A, VectorSpace B>
bool are_dual() {
    return Dual<A,B>;
};

template <typename T, size_t n>
class Rn {
    public:
    using type = Eigen::Vector<T,n>;
    static bool contains(type& x) { return true; }
};

template <typename T, size_t n>
class EuclideanSpace {
    public:
    using F = RealNumberField<T>;
    using S = F::Set;
    using V = Rn<T,n>;
    using Vec = V::type;
    using Scalar = F::Set::type;
    class Add {
        public:
        using Domain = CartesianProduct<V,V>;
        using Target = V;

        static inline Vec identity() { return Vec::Zero(); }
        static consteval bool has_identity(const Vec& x) { 
            return apply(identity(), x) == x;
        }

        static inline Vec inverse(const Vec& x) { return -x; }
        static consteval bool has_inverse(const Vec& x) { 
            return apply(inverse(x),x) == identity();
        }

        static inline Vec apply(const Vec& x, const Vec& y) { return x + y; }
        static inline Vec apply(const Domain::type& x) { return apply(x.first, x.second); }

        static consteval bool is_commutative(const Vec& x, const Vec& y) {
            return apply(x,y) == apply(y,x);
        }
        static constexpr bool is_associative(const Vec& x, const Vec& y, const Vec& z) {
            return apply(x,apply(y,z)) == apply(apply(x,y),z);
        }
    };
    class Mul {
        public:
        using Domain = CartesianProduct<S, V>;
        using Target = V;

        static inline Vec apply(const Scalar& s, const Vec& v) { return s*v; }
        static inline Vec apply(const Domain::type& x) { return apply(x.first, x.second); }
        
        static constexpr bool distributes_through(const Scalar& a, const Vec& x, const Vec& y) {
            return apply(a,Add::apply(x,y)) == Add::apply(apply(a,x),apply(a,y));
        }
        static constexpr bool distributes_through(const Scalar& a, const Scalar& b, const Vec& x) {
            return apply(F::Add::apply(a,b), x) == apply(a,x) + apply(b,x);
        }
    };
    class Contraction {
        public:
        using Domain = CartesianProduct<V,V>;
        using Target = S;
        static inline Scalar apply(const Vec& covector, const Vec& vector) { return dot(covector,vector); }
    };
};
