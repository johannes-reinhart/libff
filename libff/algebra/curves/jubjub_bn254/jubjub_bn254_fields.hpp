#ifndef JUBJUB_BN254_FIELDS_HPP_
#define JUBJUB_BN254_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t jubjub_bn254_r_bitcount = 251;
const mp_size_t jubjub_bn254_q_bitcount = 254;

const mp_size_t jubjub_bn254_r_limbs = (jubjub_bn254_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t jubjub_bn254_q_limbs = (jubjub_bn254_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<jubjub_bn254_r_limbs> jubjub_bn254_modulus_r;
extern bigint<jubjub_bn254_q_limbs> jubjub_bn254_modulus_q;

typedef Fp_model<jubjub_bn254_r_limbs, jubjub_bn254_modulus_r> jubjub_bn254_Fr;
typedef Fp_model<jubjub_bn254_q_limbs, jubjub_bn254_modulus_q> jubjub_bn254_Fq;
typedef Fp3_model<jubjub_bn254_q_limbs, jubjub_bn254_modulus_q> jubjub_bn254_Fq3;
typedef Fp6_2over3_model<jubjub_bn254_q_limbs, jubjub_bn254_modulus_q> jubjub_bn254_Fq6;
typedef jubjub_bn254_Fq6 jubjub_bn254_GT;

void init_jubjub_bn254_fields();

} // namespace libff
#endif // JUBJUB_BN254_FIELDS_HPP_
