#pragma once
#include <exception>
#include "groups.hpp"

template <class F>
concept Field = 
        AbelianGroup<typename F::Set, typename F::Add> &&
        AbelianGroup<typename F::Set, typename F::Mul> &&
requires(F::Set::type& x, F::Set::type& y, F::Set::type& z) {
    F::Mul::distributes_through<F::Add>(x,y,z);
};

template <class F>
bool is_field() {
    return Field<F>;
}

template <typename T>
class RealNumberField {
public:
    using Set = RealNumbers<T>;
    using Add = RealAddition<T>;
    using Mul = RealMultiplication<T>;
};