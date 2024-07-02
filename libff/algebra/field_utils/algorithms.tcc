/** @file
 *****************************************************************************
 Declaration of interfaces for (square-and-multiply) exponentiation and
 Tonelli-Shanks square root.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef ALGORITHMS_TCC_
#define ALGORITHMS_TCC_

#include "libff/common/utils.hpp"
#include "libff/common/profiling.hpp"

namespace libff {

using std::size_t;

template<typename FieldT, mp_size_t m>
FieldT power(const FieldT &base, const bigint<m> &exponent)
{
    FieldT result = FieldT::one();
    bool found_one = false;

    for (long i = exponent.max_bits() - 1; i >= 0; --i)
    {
        if (found_one)
        {
            result = result * result;
        }

        if (exponent.test_bit(i))
        {
            found_one = true;
            result = result * base;
        }
    }

    return result;
}

template<typename FieldT>
FieldT power(const FieldT &base, const unsigned long exponent)
{
    return power<FieldT>(base, bigint<1>(exponent));
}

template<typename FieldT>
FieldT power(const FieldT &base, const unsigned long long exponent)
{
    FieldT result = FieldT::one();

    bool found_one = false;

    for (long i = 8 * sizeof(exponent) - 1; i >= 0; --i)
    {
        if (found_one)
        {
            result = result.squared();
        }

        if (exponent & (1ull << i))
        {
            found_one = true;
            result *= base;
        }
    }

    return result;
}

template<typename FieldT>
FieldT power(const FieldT &base, const std::vector<unsigned long long> exponent)
{
    FieldT result = FieldT::one();

    bool found_one = false;

    for (unsigned long long j = 0; j < exponent.size(); j++)
    {
        unsigned long long cur_exp = exponent[j];
        for (long i = 8 * sizeof(cur_exp) - 1; i >= 0; --i)
        {
            if (found_one)
            {
                result = result.squared();
            }

            if (cur_exp & (1ull << i))
            {
                found_one = true;
                result *= base;
            }
        }
    }

    return result;
}

template<typename FieldT>
FieldT tonelli_shanks_sqrt(const FieldT &value)
{
    // A few assertions to make sure s, t, and nqr are initialized.
    assert(FieldT::s != 0);
    assert(!FieldT::t.is_even()); // Check that t is odd.
    assert(!FieldT::nqr.is_zero());

    if (value.is_zero())
    {
        return FieldT::zero();
    }

    FieldT one = FieldT::one();

    size_t v = FieldT::s;
    FieldT z = FieldT::nqr_to_t;
    FieldT w = value^FieldT::t_minus_1_over_2;
    FieldT x = value * w;
    FieldT b = x * w; // b = value^t

#if DEBUG
    // check if square with euler's criterion
    FieldT check = b;
    for (size_t i = 0; i < v-1; ++i)
    {
        check = check.squared();
    }
    assert(check == one);
#endif

    // compute square root with Tonelli--Shanks
    // (does not terminate if not a square!)

    while (b != one)
    {
        size_t m = 0;
        FieldT b2m = b;
        while (b2m != one)
        {
            /* invariant: b2m = b^(2^m) after entering this loop */
            b2m = b2m.squared();
            m += 1;
        }

        int j = v-m-1;
        w = z;
        while (j > 0)
        {
            w = w.squared();
            --j;
        } // w = z^2^(v-m-1)

        z = w.squared();
        b = b * z;
        x = x * w;
        v = m;
    }

    return x;
}

/**
 * Jacobi Symbol
 * Algorithm 2.149 from Handbook of Applied Cryptography, A. Menezes
 */
// TODO: test this
template<mp_size_t bn>
int jacobi(bigint<bn> a, bigint<bn> n){

    //Input an odd integer n ≥ 3
    assert(!(n < 3));
    assert(!n.is_even());

    //1. If a =0 then return(0)
    if (a.is_zero()){
        return 0;
    }

    //2. If a =1 then return(1)
    if (a == 1UL){
        return 1;
    }

    //3. Write a =2^e*a1, where a1 is odd.
    bigint<bn> a1 = a;
    int e = 0;
    int s = 0;
    while (a1.is_even()){
        a1 >>= 1; //a1 = a1 / 2;
        ++e;
    }

    //4. If e is even then set s←1. Otherwise set s←1 ifn ≡ 1 or 7(mod 8),
    if (e % 2 == 0 || n.as_ulong() % 8 == 1 || n.as_ulong() % 8 == 7){ // as limb is multiple of 8, we can just take the least significant limb for modulo 8
        s = 1;
    }else{ // or set s←−1 if n ≡ 3 or 5(mod 8) [Other cases do not exist, as n is odd]
        s = -1;
    }

    if (a1 == 1){ //7. If a1 =1 then return(s);
        return s;
    }
    // 7. ... otherwise return(s · JACOBI(n1,a1)).

    if (n.as_ulong() % 4 == 3 && a1.as_ulong() % 4 == 3){ //5. If n ≡ 3(mod 4) and a1 ≡ 3(mod 4) then set s←−s. [This case does not include a1=1]
        s = -s;
    }

    //6. Set n1←n mod a1. [and 7., s. above]
    return s * jacobi( n % a1, a1 );
}


template<typename FieldT>
int jacobi(const FieldT &value){
    const auto n = value.mod; // bigint type
    const auto a = value.as_bigint();
    return jacobi(a, n);
}

} // namespace libff

#endif // ALGORITHMS_TCC_
