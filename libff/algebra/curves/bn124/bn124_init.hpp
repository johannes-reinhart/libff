#ifndef BN124_INIT_HPP_
#define BN124_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/bn124/bn124_fields.hpp>

namespace libff {

// parameters for Barreto--Naehrig curve E/Fq : y^2 = x^3 + b
extern bn124_Fq bn124_coeff_b;
// parameters for twisted Barreto--Naehrig curve E'/Fq2 : y^2 = x^3 + b/xi
extern bn124_Fq2 bn124_twist;
extern bn124_Fq2 bn124_twist_coeff_b;
extern bn124_Fq bn124_twist_mul_by_b_c0;
extern bn124_Fq bn124_twist_mul_by_b_c1;
extern bn124_Fq2 bn124_twist_mul_by_q_X;
extern bn124_Fq2 bn124_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<bn124_q_limbs> bn124_ate_loop_count;
extern bool bn124_ate_is_loop_count_neg;
extern bigint<12*bn124_q_limbs> bn124_final_exponent;
extern bigint<bn124_q_limbs> bn124_final_exponent_z;
extern bool bn124_final_exponent_is_z_neg;

void init_bn124_params();

class bn124_G1;
class bn124_G2;

} // namespace libff
#endif // BN124_INIT_HPP_
