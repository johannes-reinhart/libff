#ifndef EDWARDS181_INIT_HPP_
#define EDWARDS181_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/edwards181/edwards181_fields.hpp>

namespace libff {

// parameters for Edwards curve E_{a,d}(F_q)
// parameters for twisted Edwards curve E_{a',d'}(F_q^3)
extern edwards181_Fq3 edwards181_twist;
extern edwards181_Fq3 edwards181_twist_coeff_a;
extern edwards181_Fq3 edwards181_twist_coeff_d;
extern edwards181_Fq edwards181_twist_mul_by_a_c0;
extern edwards181_Fq edwards181_twist_mul_by_a_c1;
extern edwards181_Fq edwards181_twist_mul_by_a_c2;
extern edwards181_Fq edwards181_twist_mul_by_d_c0;
extern edwards181_Fq edwards181_twist_mul_by_d_c1;
extern edwards181_Fq edwards181_twist_mul_by_d_c2;
extern edwards181_Fq edwards181_twist_mul_by_q_Y;
extern edwards181_Fq edwards181_twist_mul_by_q_Z;

// parameters for pairing
extern bigint<edwards181_q_limbs> edwards181_ate_loop_count;
extern bigint<edwards181_q_limbs> edwards181_final_exponent_last_chunk_abs_of_w0;
extern bool edwards181_final_exponent_last_chunk_is_w0_neg;
extern bigint<edwards181_q_limbs> edwards181_final_exponent_last_chunk_w1;

void init_edwards181_params();

class edwards181_G1;
class edwards181_G2;

} // libff
#endif // EDWARDS181_INIT_HPP_
