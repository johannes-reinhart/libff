#ifndef JUBJUB_ED97_FIELDS_HPP_
#define JUBJUB_ED97_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t jubjub_ed97_r_bitcount = 94;
const mp_size_t jubjub_ed97_q_bitcount = 97;

const mp_size_t jubjub_ed97_r_limbs = (jubjub_ed97_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t jubjub_ed97_q_limbs = (jubjub_ed97_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<jubjub_ed97_r_limbs> jubjub_ed97_modulus_r;
extern bigint<jubjub_ed97_q_limbs> jubjub_ed97_modulus_q;

typedef Fp_model<jubjub_ed97_r_limbs, jubjub_ed97_modulus_r> jubjub_ed97_Fr;
typedef Fp_model<jubjub_ed97_q_limbs, jubjub_ed97_modulus_q> jubjub_ed97_Fq;
typedef Fp3_model<jubjub_ed97_q_limbs, jubjub_ed97_modulus_q> jubjub_ed97_Fq3;
typedef Fp6_2over3_model<jubjub_ed97_q_limbs, jubjub_ed97_modulus_q> jubjub_ed97_Fq6;
typedef jubjub_ed97_Fq6 jubjub_ed97_GT;

void init_jubjub_ed97_fields();

} // namespace libff
#endif // JUBJUB_ED97_FIELDS_HPP_
