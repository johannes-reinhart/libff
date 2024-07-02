#ifndef EDWARDS97_FIELDS_HPP_
#define EDWARDS97_FIELDS_HPP_
#include <libff/algebra/fields/prime_base/fp.hpp>
#include <libff/algebra/fields/prime_extension/fp3.hpp>
#include <libff/algebra/fields/prime_extension/fp6_2over3.hpp>

namespace libff {

const mp_size_t edwards97_r_bitcount = 97;
const mp_size_t edwards97_q_bitcount = 99;

const mp_size_t edwards97_r_limbs = (edwards97_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t edwards97_q_limbs = (edwards97_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<edwards97_r_limbs> edwards97_modulus_r;
extern bigint<edwards97_q_limbs> edwards97_modulus_q;

typedef Fp_model<edwards97_r_limbs, edwards97_modulus_r> edwards97_Fr;
typedef Fp_model<edwards97_q_limbs, edwards97_modulus_q> edwards97_Fq;
typedef Fp3_model<edwards97_q_limbs, edwards97_modulus_q> edwards97_Fq3;
typedef Fp6_2over3_model<edwards97_q_limbs, edwards97_modulus_q> edwards97_Fq6;
typedef edwards97_Fq6 edwards97_GT;

void init_edwards97_fields();

} // namespace libff
#endif // EDWARDS97_FIELDS_HPP_
