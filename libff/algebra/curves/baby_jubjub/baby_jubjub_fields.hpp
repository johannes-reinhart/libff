#ifndef BABY_JUBJUB_FIELDS_HPP_
#define BABY_JUBJUB_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t baby_jubjub_r_bitcount = 251;
const mp_size_t baby_jubjub_q_bitcount = 254;

const mp_size_t baby_jubjub_r_limbs = (baby_jubjub_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t baby_jubjub_q_limbs = (baby_jubjub_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<baby_jubjub_r_limbs> baby_jubjub_modulus_r;
extern bigint<baby_jubjub_q_limbs> baby_jubjub_modulus_q;

typedef Fp_model<baby_jubjub_r_limbs, baby_jubjub_modulus_r> baby_jubjub_Fr;
typedef Fp_model<baby_jubjub_q_limbs, baby_jubjub_modulus_q> baby_jubjub_Fq;
typedef Fp3_model<baby_jubjub_q_limbs, baby_jubjub_modulus_q> baby_jubjub_Fq3;
typedef Fp6_2over3_model<baby_jubjub_q_limbs, baby_jubjub_modulus_q> baby_jubjub_Fq6;
typedef baby_jubjub_Fq6 baby_jubjub_GT;

void init_baby_jubjub_fields();

} // namespace libff
#endif // BABY_JUBJUB_FIELDS_HPP_
