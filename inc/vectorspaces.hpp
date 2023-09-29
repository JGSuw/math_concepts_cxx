#pragma once
#include "fields.hpp"
#include <Eigen/Dense>


template <class VecSp>
concept VectorSpace = 
        Field<typename VecSp::F> &&
        AbelianGroup<typename VecSp::V, typename VecSp::Add> &&
        MapsTo<
            typename VecSp::Mul,
            CartesianProduct<typename VecSp::F::Set, typename VecSp::V>, 
            typename VecSp::V> &&
requires (VecSp::F::Set::type& a, VecSp::F::Set::type& b, VecSp::V::type& x, VecSp::V::type& y) {
    VecSp::Mul::unity_scalar_fixes_identity(x);
    VecSp::Mul::distributes_through_vec_add(a,x,y);
    VecSp::Mul::distributes_through_scalar_add(a,b,x);
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
        static void has_identity(const Vec& x) { 
            static_assert(apply(identity(), x) == x);
        }

        static inline Vec inverse(const Vec& x) { return -x; }
        static void has_inverse(const Vec& x) { 
            static_assert(apply(inverse(x),x) == identity());
        }

        static inline Vec apply(const Vec& x, const Vec& y) { return x + y; }
        static inline Vec apply(const Domain::type& x) { return apply(x.first, x.second); }

        static void is_commutative(const Vec& x, const Vec& y) {
            static_assert(apply(x,y) == apply(y,x));
        }
        static void is_associative(const Vec& x, const Vec& y, const Vec& z) {
            static_assert(apply(x,apply(y,z)) == apply(apply(x,y),z));
        }
    };
    class Mul {
        public:
        using Domain = CartesianProduct<S, V>;
        using Target = V;

        static inline Vec apply(const Scalar& s, const Vec& v) { return s*v; }
        static inline Vec apply(const Domain::type& x) { return apply(x.first, x.second); }

        static void unity_scalar_fixes_identity(const Vec& x) {
            static_assert(apply(F::Mul::identity(), x) == x)
        }
  
        static void distributes_through_vec_add(const Scalar& a, const Vec& x, const Vec& y) {
            static_assert(apply(a,Add::apply(x,y)) == Add::apply(apply(a,x),apply(a,y)));
        }

        static void distributes_through_scalar_add(const Scalar& a, const Scalar& b, const Vec& x) {
            static_assert(apply(F::Add::apply(a,b), x) == apply(a,x) + apply(b,x));
        }

    };
    class Contraction {
        public:
        using Domain = CartesianProduct<V,V>;
        using Target = S;
        static inline Scalar apply(const Vec& covector, const Vec& vector) { return dot(covector,vector); }
    };
};
