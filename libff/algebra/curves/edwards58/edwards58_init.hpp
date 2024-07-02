#ifndef EDWARDS58_INIT_HPP_
#define EDWARDS58_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/edwards58/edwards58_fields.hpp>

namespace libff {

// parameters for Edwards curve E_{a,d}(F_q)
// parameters for twisted Edwards curve E_{a',d'}(F_q^3)
extern edwards58_Fq3 edwards58_twist;
extern edwards58_Fq3 edwards58_twist_coeff_a;
extern edwards58_Fq3 edwards58_twist_coeff_d;
extern edwards58_Fq edwards58_twist_mul_by_a_c0;
extern edwards58_Fq edwards58_twist_mul_by_a_c1;
extern edwards58_Fq edwards58_twist_mul_by_a_c2;
extern edwards58_Fq edwards58_twist_mul_by_d_c0;
extern edwards58_Fq edwards58_twist_mul_by_d_c1;
extern edwards58_Fq edwards58_twist_mul_by_d_c2;
extern edwards58_Fq edwards58_twist_mul_by_q_Y;
extern edwards58_Fq edwards58_twist_mul_by_q_Z;

// parameters for pairing
extern bigint<edwards58_q_limbs> edwards58_ate_loop_count;
extern bigint<edwards58_q_limbs> edwards58_final_exponent_last_chunk_abs_of_w0;
extern bool edwards58_final_exponent_last_chunk_is_w0_neg;
extern bigint<edwards58_q_limbs> edwards58_final_exponent_last_chunk_w1;

void init_edwards58_params();

class edwards58_G1;
class edwards58_G2;

} // libff
#endif // EDWARDS58_INIT_HPP_
