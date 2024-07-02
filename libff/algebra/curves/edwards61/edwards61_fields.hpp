#ifndef EDWARDS61_FIELDS_HPP_
#define EDWARDS61_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t edwards61_r_bitcount = 61;
const mp_size_t edwards61_q_bitcount = 63;

const mp_size_t edwards61_r_limbs = (edwards61_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t edwards61_q_limbs = (edwards61_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<edwards61_r_limbs> edwards61_modulus_r;
extern bigint<edwards61_q_limbs> edwards61_modulus_q;

typedef Fp_model<edwards61_r_limbs, edwards61_modulus_r> edwards61_Fr;
typedef Fp_model<edwards61_q_limbs, edwards61_modulus_q> edwards61_Fq;
typedef Fp3_model<edwards61_q_limbs, edwards61_modulus_q> edwards61_Fq3;
typedef Fp6_2over3_model<edwards61_q_limbs, edwards61_modulus_q> edwards61_Fq6;
typedef edwards61_Fq6 edwards61_GT;

void init_edwards61_fields();

} // namespace libff
#endif // EDWARDS61_FIELDS_HPP_
