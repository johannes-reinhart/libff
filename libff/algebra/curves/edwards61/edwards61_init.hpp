#ifndef EDWARDS61_INIT_HPP_
#define EDWARDS61_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/edwards61/edwards61_fields.hpp>

namespace libff {

// parameters for Edwards curve E_{a,d}(F_q)
// parameters for twisted Edwards curve E_{a',d'}(F_q^3)
extern edwards61_Fq3 edwards61_twist;
extern edwards61_Fq3 edwards61_twist_coeff_a;
extern edwards61_Fq3 edwards61_twist_coeff_d;
extern edwards61_Fq edwards61_twist_mul_by_a_c0;
extern edwards61_Fq edwards61_twist_mul_by_a_c1;
extern edwards61_Fq edwards61_twist_mul_by_a_c2;
extern edwards61_Fq edwards61_twist_mul_by_d_c0;
extern edwards61_Fq edwards61_twist_mul_by_d_c1;
extern edwards61_Fq edwards61_twist_mul_by_d_c2;
extern edwards61_Fq edwards61_twist_mul_by_q_Y;
extern edwards61_Fq edwards61_twist_mul_by_q_Z;

// parameters for pairing
extern bigint<edwards61_q_limbs> edwards61_ate_loop_count;
extern bigint<edwards61_q_limbs> edwards61_final_exponent_last_chunk_abs_of_w0;
extern bool edwards61_final_exponent_last_chunk_is_w0_neg;
extern bigint<edwards61_q_limbs> edwards61_final_exponent_last_chunk_w1;

void init_edwards61_params();

class edwards61_G1;
class edwards61_G2;

} // libff
#endif // EDWARDS61_INIT_HPP_
