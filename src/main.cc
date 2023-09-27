#include <iostream>
#include <exception>
// #include "fields.hpp"
#include "vectorspaces.hpp"

using R = RealNumbers<double>;
using Add = RealAddition<double>;
using Mul = RealMultiplication<double>;
using Reals = RealNumberField<double>;
using V = Rn<double,3>;
using VecSp = EuclideanSpace<double,3>;

int main(void) {
    is_abelian_group<R,Add>();
    is_abelian_group<R,Mul>();
    is_field<Reals>();
    is_pointset<Reals::Set>();
    is_pointset<V>();
    is_vectorspace<VecSp>();
    // are_dual<VecSp,VecSp>();
    return 0;
}