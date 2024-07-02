#ifndef JUBJUB_BN183_FIELDS_HPP_
#define JUBJUB_BN183_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t jubjub_bn183_r_bitcount = 180;
const mp_size_t jubjub_bn183_q_bitcount = 183;

const mp_size_t jubjub_bn183_r_limbs = (jubjub_bn183_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t jubjub_bn183_q_limbs = (jubjub_bn183_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<jubjub_bn183_r_limbs> jubjub_bn183_modulus_r;
extern bigint<jubjub_bn183_q_limbs> jubjub_bn183_modulus_q;

typedef Fp_model<jubjub_bn183_r_limbs, jubjub_bn183_modulus_r> jubjub_bn183_Fr;
typedef Fp_model<jubjub_bn183_q_limbs, jubjub_bn183_modulus_q> jubjub_bn183_Fq;
typedef Fp3_model<jubjub_bn183_q_limbs, jubjub_bn183_modulus_q> jubjub_bn183_Fq3;
typedef Fp6_2over3_model<jubjub_bn183_q_limbs, jubjub_bn183_modulus_q> jubjub_bn183_Fq6;
typedef jubjub_bn183_Fq6 jubjub_bn183_GT;

void init_jubjub_bn183_fields();

} // namespace libff
#endif // JUBJUB_BN183_FIELDS_HPP_
