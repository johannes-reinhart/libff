#ifndef EDWARDS181_FIELDS_HPP_
#define EDWARDS181_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t edwards181_r_bitcount = 181;
const mp_size_t edwards181_q_bitcount = 183;

const mp_size_t edwards181_r_limbs = (edwards181_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t edwards181_q_limbs = (edwards181_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<edwards181_r_limbs> edwards181_modulus_r;
extern bigint<edwards181_q_limbs> edwards181_modulus_q;

typedef Fp_model<edwards181_r_limbs, edwards181_modulus_r> edwards181_Fr;
typedef Fp_model<edwards181_q_limbs, edwards181_modulus_q> edwards181_Fq;
typedef Fp3_model<edwards181_q_limbs, edwards181_modulus_q> edwards181_Fq3;
typedef Fp6_2over3_model<edwards181_q_limbs, edwards181_modulus_q> edwards181_Fq6;
typedef edwards181_Fq6 edwards181_GT;

void init_edwards181_fields();

} // namespace libff
#endif // EDWARDS181_FIELDS_HPP_
