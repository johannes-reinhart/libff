/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN254_FIELDS_HPP_
#define BN254_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp12_2over3over2.hpp>
#include <libff/algebra/fields/prime_extension/fp2.hpp>
#include <libff/algebra/fields/prime_extension/fp6_3over2.hpp>

namespace libff {

const mp_size_t bn254_r_bitcount = 254;
const mp_size_t bn254_q_bitcount = 254;

const mp_size_t bn254_r_limbs = (bn254_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t bn254_q_limbs = (bn254_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<bn254_r_limbs> bn254_modulus_r;
extern bigint<bn254_q_limbs> bn254_modulus_q;

typedef Fp_model<bn254_r_limbs, bn254_modulus_r> bn254_Fr;
typedef Fp_model<bn254_q_limbs, bn254_modulus_q> bn254_Fq;
typedef Fp2_model<bn254_q_limbs, bn254_modulus_q> bn254_Fq2;
typedef Fp6_3over2_model<bn254_q_limbs, bn254_modulus_q> bn254_Fq6;
typedef Fp12_2over3over2_model<bn254_q_limbs, bn254_modulus_q> bn254_Fq12;
typedef bn254_Fq12 bn254_GT;

void init_bn254_fields();

} // namespace libff

#endif // BN254_FIELDS_HPP_
