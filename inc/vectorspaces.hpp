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
requires (A::Add& addition) {
    A::Mul::distributes_through(addition);
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
        static consteval Vec identity() { return Vec::Zero(); }
        static consteval Vec inverse(const Vec& x) { return -x; }
        static constexpr bool is_commutative = true;
        static constexpr bool is_associative = true;
        static consteval Vec apply(const Vec& x, const Vec& y) { return x + y; }
        static consteval Vec apply(const Domain::type& x) { return apply(x.first, x.second); }
    };
    class Mul {
        public:
        using Domain = CartesianProduct<S, V>;
        using Target = V;
        static consteval Vec apply(const Scalar& s, const Vec& v) { return s*v; }
        static consteval Vec apply(const Domain::type& x) { return apply(x.first, x.second); }
        static constexpr bool distributes_through(EuclideanSpace<T,n>::Add& add) {return true;}
    };
    class Contraction {
        public:
        using Domain = CartesianProduct<V,V>;
        using Target = S;
        static consteval Scalar apply(const Vec& covector, const Vec& vector) { return dot(covector,vector); }
    };
};
