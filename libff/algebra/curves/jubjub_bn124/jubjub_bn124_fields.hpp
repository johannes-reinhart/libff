#ifndef JUBJUB_BN124_FIELDS_HPP_
#define JUBJUB_BN124_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t jubjub_bn124_r_bitcount = 121;
const mp_size_t jubjub_bn124_q_bitcount = 124;

const mp_size_t jubjub_bn124_r_limbs = (jubjub_bn124_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t jubjub_bn124_q_limbs = (jubjub_bn124_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<jubjub_bn124_r_limbs> jubjub_bn124_modulus_r;
extern bigint<jubjub_bn124_q_limbs> jubjub_bn124_modulus_q;

typedef Fp_model<jubjub_bn124_r_limbs, jubjub_bn124_modulus_r> jubjub_bn124_Fr;
typedef Fp_model<jubjub_bn124_q_limbs, jubjub_bn124_modulus_q> jubjub_bn124_Fq;
typedef Fp3_model<jubjub_bn124_q_limbs, jubjub_bn124_modulus_q> jubjub_bn124_Fq3;
typedef Fp6_2over3_model<jubjub_bn124_q_limbs, jubjub_bn124_modulus_q> jubjub_bn124_Fq6;
typedef jubjub_bn124_Fq6 jubjub_bn124_GT;

void init_jubjub_bn124_fields();

} // namespace libff
#endif // JUBJUB_BN124_FIELDS_HPP_
