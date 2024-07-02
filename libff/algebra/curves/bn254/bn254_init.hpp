#ifndef BN254_INIT_HPP_
#define BN254_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/bn254/bn254_fields.hpp>

namespace libff {

// parameters for Barreto--Naehrig curve E/Fq : y^2 = x^3 + b
extern bn254_Fq bn254_coeff_b;
// parameters for twisted Barreto--Naehrig curve E'/Fq2 : y^2 = x^3 + b/xi
extern bn254_Fq2 bn254_twist;
extern bn254_Fq2 bn254_twist_coeff_b;
extern bn254_Fq bn254_twist_mul_by_b_c0;
extern bn254_Fq bn254_twist_mul_by_b_c1;
extern bn254_Fq2 bn254_twist_mul_by_q_X;
extern bn254_Fq2 bn254_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<bn254_q_limbs> bn254_ate_loop_count;
extern bool bn254_ate_is_loop_count_neg;
extern bigint<12*bn254_q_limbs> bn254_final_exponent;
extern bigint<bn254_q_limbs> bn254_final_exponent_z;
extern bool bn254_final_exponent_is_z_neg;

void init_bn254_params();

class bn254_G1;
class bn254_G2;

} // namespace libff
#endif // BN254_INIT_HPP_
