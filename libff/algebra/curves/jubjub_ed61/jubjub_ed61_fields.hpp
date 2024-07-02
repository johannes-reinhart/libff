#ifndef JUBJUB_ED61_FIELDS_HPP_
#define JUBJUB_ED61_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t jubjub_ed61_r_bitcount = 58;
const mp_size_t jubjub_ed61_q_bitcount = 61;

const mp_size_t jubjub_ed61_r_limbs = (jubjub_ed61_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t jubjub_ed61_q_limbs = (jubjub_ed61_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<jubjub_ed61_r_limbs> jubjub_ed61_modulus_r;
extern bigint<jubjub_ed61_q_limbs> jubjub_ed61_modulus_q;

typedef Fp_model<jubjub_ed61_r_limbs, jubjub_ed61_modulus_r> jubjub_ed61_Fr;
typedef Fp_model<jubjub_ed61_q_limbs, jubjub_ed61_modulus_q> jubjub_ed61_Fq;
typedef Fp3_model<jubjub_ed61_q_limbs, jubjub_ed61_modulus_q> jubjub_ed61_Fq3;
typedef Fp6_2over3_model<jubjub_ed61_q_limbs, jubjub_ed61_modulus_q> jubjub_ed61_Fq6;
typedef jubjub_ed61_Fq6 jubjub_ed61_GT;

void init_jubjub_ed61_fields();

} // namespace libff
#endif // JUBJUB_ED61_FIELDS_HPP_
