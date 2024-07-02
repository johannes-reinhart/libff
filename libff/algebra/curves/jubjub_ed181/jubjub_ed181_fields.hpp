#ifndef JUBJUB_ED181_FIELDS_HPP_
#define JUBJUB_ED181_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t jubjub_ed181_r_bitcount = 178;
const mp_size_t jubjub_ed181_q_bitcount = 181;

const mp_size_t jubjub_ed181_r_limbs = (jubjub_ed181_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t jubjub_ed181_q_limbs = (jubjub_ed181_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<jubjub_ed181_r_limbs> jubjub_ed181_modulus_r;
extern bigint<jubjub_ed181_q_limbs> jubjub_ed181_modulus_q;

typedef Fp_model<jubjub_ed181_r_limbs, jubjub_ed181_modulus_r> jubjub_ed181_Fr;
typedef Fp_model<jubjub_ed181_q_limbs, jubjub_ed181_modulus_q> jubjub_ed181_Fq;
typedef Fp3_model<jubjub_ed181_q_limbs, jubjub_ed181_modulus_q> jubjub_ed181_Fq3;
typedef Fp6_2over3_model<jubjub_ed181_q_limbs, jubjub_ed181_modulus_q> jubjub_ed181_Fq6;
typedef jubjub_ed181_Fq6 jubjub_ed181_GT;

void init_jubjub_ed181_fields();

} // namespace libff
#endif // JUBJUB_ED181_FIELDS_HPP_
