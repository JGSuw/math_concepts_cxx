#pragma once
#include <Eigen/Dense>
#include "operators.hpp"
using namespace Eigen;

template <class T>
concept PointSet = requires(T::type& x, bool& b) {
    b = T::contains(x);
};

template <class T>
bool is_pointset() {
    return PointSet<T>;
}

template <PointSet X, PointSet Y>
class CartesianProduct{
    public:
    using type = std::pair<typename X::type, typename Y::type>;
    static consteval bool contains(X::type& x, Y::type& y) {
        return X::contains(x) && Y::contains(y);
    }
    static consteval bool contains(type& pair) {
        return contains(pair.first, pair.second);
    }
};

template <typename T>
class NaturalNumbers {
    public:
    using type = T;
    static constexpr T zero = static_cast<T>(0);
    static constexpr T one = static_cast<T>(1);
    static consteval bool contains(T& x) {
        return x >= zero && x % one == zero;
    };
};

template <typename T>
class Integers {
    public:
    using type = T;
    static constexpr T zero = NaturalNumbers<T>::zero;
    static constexpr T one = NaturalNumbers<T>::one;
    static consteval bool contains(T& x) {
        return x % one == zero;
    };
};

template <typename T>
class RationalNumbers {
    public:
    using type = std::pair<T,T>;
    static consteval bool contains(T& x)  {
        return (Integers<T>::contains(x.first) &&
                Integers<T>::contains(x.second) &&
                x.second != Integers<T>::zero);
    }
};

template <typename T>
class RealNumbers {
    public:
    using type = T;
    static constexpr T zero = NaturalNumbers<T>::zero;
    static constexpr T one = NaturalNumbers<T>::one;
    static consteval bool contains(T& x) { 
        return true;
    }
};