/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN183_FIELDS_HPP_
#define BN183_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp12_2over3over2.hpp>
#include <libff/algebra/fields/prime_extension/fp2.hpp>
#include <libff/algebra/fields/prime_extension/fp6_3over2.hpp>

namespace libff {

const mp_size_t bn183_r_bitcount = 183;
const mp_size_t bn183_q_bitcount = 183;

const mp_size_t bn183_r_limbs = (bn183_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t bn183_q_limbs = (bn183_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<bn183_r_limbs> bn183_modulus_r;
extern bigint<bn183_q_limbs> bn183_modulus_q;

typedef Fp_model<bn183_r_limbs, bn183_modulus_r> bn183_Fr;
typedef Fp_model<bn183_q_limbs, bn183_modulus_q> bn183_Fq;
typedef Fp2_model<bn183_q_limbs, bn183_modulus_q> bn183_Fq2;
typedef Fp6_3over2_model<bn183_q_limbs, bn183_modulus_q> bn183_Fq6;
typedef Fp12_2over3over2_model<bn183_q_limbs, bn183_modulus_q> bn183_Fq12;
typedef bn183_Fq12 bn183_GT;

void init_bn183_fields();

} // namespace libff

#endif // BN183_FIELDS_HPP_
