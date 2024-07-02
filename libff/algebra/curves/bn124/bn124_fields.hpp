/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN124_FIELDS_HPP_
#define BN124_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp12_2over3over2.hpp>
#include <libff/algebra/fields/prime_extension/fp2.hpp>
#include <libff/algebra/fields/prime_extension/fp6_3over2.hpp>

namespace libff {

const mp_size_t bn124_r_bitcount = 124;
const mp_size_t bn124_q_bitcount = 124;

const mp_size_t bn124_r_limbs = (bn124_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t bn124_q_limbs = (bn124_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<bn124_r_limbs> bn124_modulus_r;
extern bigint<bn124_q_limbs> bn124_modulus_q;

typedef Fp_model<bn124_r_limbs, bn124_modulus_r> bn124_Fr;
typedef Fp_model<bn124_q_limbs, bn124_modulus_q> bn124_Fq;
typedef Fp2_model<bn124_q_limbs, bn124_modulus_q> bn124_Fq2;
typedef Fp6_3over2_model<bn124_q_limbs, bn124_modulus_q> bn124_Fq6;
typedef Fp12_2over3over2_model<bn124_q_limbs, bn124_modulus_q> bn124_Fq12;
typedef bn124_Fq12 bn124_GT;

void init_bn124_fields();

} // namespace libff

#endif // BN124_FIELDS_HPP_
