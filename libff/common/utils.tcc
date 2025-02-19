/** @file
 *****************************************************************************
 Implementation of templatized utility functions.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef UTILS_TCC_
#define UTILS_TCC_

namespace libff {

using std::size_t;

template<typename CurveT>
size_t curve_size_in_bits(const std::vector<CurveT> &v)
{
    return v.size() * CurveT::size_in_bits();
}

template<typename T>
T random_element_non_zero_one()
{
    T x = T::random_element();
    while (x.is_zero() || x == T::one())
        x = T::random_element();
    return x;
}

template<typename T>
T random_element_non_zero()
{
    T x = T::random_element();
    while (x.is_zero())
        x = T::random_element();
    return x;
}

template<typename T>
T random_element_exclude(T y)
{
    T x = T::random_element();
    while (x == y)
        x = T::random_element();
    return x;
}

template<typename Fp_to, typename Fp_from>
Fp_to convert_field(Fp_from a)
{
    assert(Fp_to::mod == Fp_from::mod); // Poseidon Parameters must mach curve parameters
    Fp_to b;
    b.mont_repr = a.mont_repr;
    return b;
}

} // namespace libff

#endif // UTILS_TCC_
