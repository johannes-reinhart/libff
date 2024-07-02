#ifndef EDWARDS58_FIELDS_HPP_
#define EDWARDS58_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t edwards58_r_bitcount = 58;
const mp_size_t edwards58_q_bitcount = 60;

const mp_size_t edwards58_r_limbs = (edwards58_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t edwards58_q_limbs = (edwards58_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<edwards58_r_limbs> edwards58_modulus_r;
extern bigint<edwards58_q_limbs> edwards58_modulus_q;

typedef Fp_model<edwards58_r_limbs, edwards58_modulus_r> edwards58_Fr;
typedef Fp_model<edwards58_q_limbs, edwards58_modulus_q> edwards58_Fq;
typedef Fp3_model<edwards58_q_limbs, edwards58_modulus_q> edwards58_Fq3;
typedef Fp6_2over3_model<edwards58_q_limbs, edwards58_modulus_q> edwards58_Fq6;
typedef edwards58_Fq6 edwards58_GT;

void init_edwards58_fields();

} // namespace libff
#endif // EDWARDS58_FIELDS_HPP_
