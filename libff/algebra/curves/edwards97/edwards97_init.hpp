#ifndef EDWARDS97_INIT_HPP_
#define EDWARDS97_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/edwards97/edwards97_fields.hpp>

namespace libff {

// parameters for Edwards curve E_{a,d}(F_q)
// parameters for twisted Edwards curve E_{a',d'}(F_q^3)
extern edwards97_Fq3 edwards97_twist;
extern edwards97_Fq3 edwards97_twist_coeff_a;
extern edwards97_Fq3 edwards97_twist_coeff_d;
extern edwards97_Fq edwards97_twist_mul_by_a_c0;
extern edwards97_Fq edwards97_twist_mul_by_a_c1;
extern edwards97_Fq edwards97_twist_mul_by_a_c2;
extern edwards97_Fq edwards97_twist_mul_by_d_c0;
extern edwards97_Fq edwards97_twist_mul_by_d_c1;
extern edwards97_Fq edwards97_twist_mul_by_d_c2;
extern edwards97_Fq edwards97_twist_mul_by_q_Y;
extern edwards97_Fq edwards97_twist_mul_by_q_Z;

// parameters for pairing
extern bigint<edwards97_q_limbs> edwards97_ate_loop_count;
extern bigint<edwards97_q_limbs> edwards97_final_exponent_last_chunk_abs_of_w0;
extern bool edwards97_final_exponent_last_chunk_is_w0_neg;
extern bigint<edwards97_q_limbs> edwards97_final_exponent_last_chunk_w1;

void init_edwards97_params();

class edwards97_G1;
class edwards97_G2;

} // libff
#endif // EDWARDS97_INIT_HPP_
