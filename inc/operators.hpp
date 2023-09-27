#pragma once
#include "point_set.hpp"
template <class F>
concept Function = 
        PointSet<typename F::Domain> && 
        PointSet<typename F::Target> && 
requires(F::Domain::type x, F::Target::type& y) {
    y = F::apply(x);
};

template <class F, class X, class Y>
concept MapsTo = Function<F> && 
        std::same_as<typename F::Domain, X> &&
        std::same_as<typename F::Target, Y>;

template <class F>
concept BinaryOperator = 
            Function<F> &&
requires (F::Domain::type& x) {
    x.first;
    x.second;
};


template <class F>
concept Commutative = BinaryOperator<F> && F::is_commutative;

template <class F>
concept ClosedBinaryOperator = BinaryOperator<F> &&
        std::same_as<typename F::Domain::type::first_type, typename F::Domain::type::second_type> &&
        std::same_as<typename F::Domain::type::first_type, typename F::Target::type>;

template <class F>
concept Associative = ClosedBinaryOperator<F> && F::is_associative;

template <class F>
concept HasIdentityElement = ClosedBinaryOperator<F> && 
requires (F::Target::type& x) {
    x = F::identity();
};

template <class F>
concept HasInverseElement = ClosedBinaryOperator<F> && 
requires (F::Target::type& x) {
    x = F::inverse(x);
};