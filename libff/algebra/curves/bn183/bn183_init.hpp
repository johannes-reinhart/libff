#ifndef BN183_INIT_HPP_
#define BN183_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/bn183/bn183_fields.hpp>

namespace libff {

// parameters for Barreto--Naehrig curve E/Fq : y^2 = x^3 + b
extern bn183_Fq bn183_coeff_b;
// parameters for twisted Barreto--Naehrig curve E'/Fq2 : y^2 = x^3 + b/xi
extern bn183_Fq2 bn183_twist;
extern bn183_Fq2 bn183_twist_coeff_b;
extern bn183_Fq bn183_twist_mul_by_b_c0;
extern bn183_Fq bn183_twist_mul_by_b_c1;
extern bn183_Fq2 bn183_twist_mul_by_q_X;
extern bn183_Fq2 bn183_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<bn183_q_limbs> bn183_ate_loop_count;
extern bool bn183_ate_is_loop_count_neg;
extern bigint<12*bn183_q_limbs> bn183_final_exponent;
extern bigint<bn183_q_limbs> bn183_final_exponent_z;
extern bool bn183_final_exponent_is_z_neg;

void init_bn183_params();

class bn183_G1;
class bn183_G2;

} // namespace libff
#endif // BN183_INIT_HPP_
