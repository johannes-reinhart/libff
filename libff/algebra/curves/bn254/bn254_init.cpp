#include <libff/algebra/curves/bn254/bn254_g1.hpp>
#include <libff/algebra/curves/bn254/bn254_g2.hpp>
#include <libff/algebra/curves/bn254/bn254_init.hpp>

namespace libff {

bn254_Fq bn254_coeff_b;
bn254_Fq2 bn254_twist;
bn254_Fq2 bn254_twist_coeff_b;
bn254_Fq bn254_twist_mul_by_b_c0;
bn254_Fq bn254_twist_mul_by_b_c1;
bn254_Fq2 bn254_twist_mul_by_q_X;
bn254_Fq2 bn254_twist_mul_by_q_Y;

bigint<bn254_q_limbs> bn254_ate_loop_count;
bool bn254_ate_is_loop_count_neg;
bigint<bn254_q_limbs> bn254_final_exponent_z;
bool bn254_final_exponent_is_z_neg;

void init_bn254_params()
{
	init_bn254_fields();

    /* choice of short Weierstrass curve and its twist */

	bn254_coeff_b = bn254_Fq("3");
	bn254_twist = bn254_Fq2(bn254_Fq("9"), bn254_Fq("1"));
	bn254_twist_coeff_b = bn254_coeff_b * bn254_twist.inverse();
	bn254_twist_mul_by_b_c0 = bn254_coeff_b * bn254_Fq2::non_residue;
	bn254_twist_mul_by_b_c1 = bn254_coeff_b * bn254_Fq2::non_residue;
	bn254_twist_mul_by_q_X = bn254_Fq2(bn254_Fq("21575463638280843010398324269430826099269044274347216827212613867836435027261"),bn254_Fq("10307601595873709700152284273816112264069230130616436755625194854815875713954"));
	bn254_twist_mul_by_q_Y = bn254_Fq2(bn254_Fq("2821565182194536844548159561693502659359617185244120367078079554186484126554"),bn254_Fq("3505843767911556378687030309984248845540243509899259641013678093033130930403"));
	
    /* choice of group G1 */

	bn254_G1::G1_zero = bn254_G1(bn254_Fq::zero(),   
                                                   bn254_Fq::one(),
                                                   bn254_Fq::zero());
	bn254_G1::G1_one = bn254_G1(bn254_Fq("1"),
                                                  bn254_Fq("2"),
                                                  bn254_Fq::one());
    bn254_G1::initialized = true;

    // Cofactor
    bn254_G1::h = bigint<bn254_G1::h_limbs>("1");

    // TODO
	bn254_G1::wnaf_window_table.resize(0);
    bn254_G1::wnaf_window_table.push_back(11);
    bn254_G1::wnaf_window_table.push_back(24);
    bn254_G1::wnaf_window_table.push_back(60);
    bn254_G1::wnaf_window_table.push_back(127);

    bn254_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.99]
    bn254_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.99, 10.99]
    bn254_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.99, 32.29]
    bn254_G1::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [32.29, 55.23]
    bn254_G1::fixed_base_exp_window_table.push_back(32);
    // window 5 is unbeaten in [55.23, 162.03]
    bn254_G1::fixed_base_exp_window_table.push_back(55);
    // window 6 is unbeaten in [162.03, 360.15]
    bn254_G1::fixed_base_exp_window_table.push_back(162);
    // window 7 is unbeaten in [360.15, 815.44]
    bn254_G1::fixed_base_exp_window_table.push_back(360);
    // window 8 is unbeaten in [815.44, 2373.07]
    bn254_G1::fixed_base_exp_window_table.push_back(815);
    // window 9 is unbeaten in [2373.07, 6977.75]
    bn254_G1::fixed_base_exp_window_table.push_back(2373);
    // window 10 is unbeaten in [6977.75, 7122.23]
    bn254_G1::fixed_base_exp_window_table.push_back(6978);
    // window 11 is unbeaten in [7122.23, 57818.46]
    bn254_G1::fixed_base_exp_window_table.push_back(7122);
    // window 12 is never the best
    bn254_G1::fixed_base_exp_window_table.push_back(0);
    // window 13 is unbeaten in [57818.46, 169679.14]
    bn254_G1::fixed_base_exp_window_table.push_back(57818);
    // window 14 is never the best
    bn254_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [169679.14, 439758.91]
    bn254_G1::fixed_base_exp_window_table.push_back(169679);
    // window 16 is unbeaten in [439758.91, 936073.41]
    bn254_G1::fixed_base_exp_window_table.push_back(439759);
    // window 17 is unbeaten in [936073.41, 4666554.74]
    bn254_G1::fixed_base_exp_window_table.push_back(936073);
    // window 18 is never the best
    bn254_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4666554.74, 7580404.42]
    bn254_G1::fixed_base_exp_window_table.push_back(4666555);
    // window 20 is unbeaten in [7580404.42, 34552892.20]
    bn254_G1::fixed_base_exp_window_table.push_back(7580404);
    // window 21 is never the best
    bn254_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [34552892.20, inf]
    bn254_G1::fixed_base_exp_window_table.push_back(34552892);

    /* choice of group G2 */

	bn254_G2::G2_zero = bn254_G2(bn254_Fq2::zero(),
                                 bn254_Fq2::one(),
                                 bn254_Fq2::zero());
	bn254_G2::G2_one = bn254_G2(bn254_Fq2(bn254_Fq("5212200097248580463769379349309309787173734510675045144092462538570320376832"),
                                                                     bn254_Fq("8347398511459179896766693881465182238445952000430972078774508851968137266123")),
                                                  bn254_Fq2(bn254_Fq("17306101117694095396656336169793795394209298809409544012716165555251634886935"),
                                                                     bn254_Fq("4680157025344042565063060947083364758580753209003414360718145377043375940606")),
                                                  bn254_Fq2::one());

    // Cofactor
    bn254_G2::h = bigint<bn254_G2::h_limbs>("21888242871839275222246405745257275088844257914179612981679871602714643921549");

    // TODO
	bn254_G2::wnaf_window_table.resize(0);
    bn254_G2::wnaf_window_table.push_back(5);
    bn254_G2::wnaf_window_table.push_back(15);
    bn254_G2::wnaf_window_table.push_back(39);
    bn254_G2::wnaf_window_table.push_back(109);

    bn254_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.10]
    bn254_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.10, 10.43]
    bn254_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.43, 25.28]
    bn254_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.28, 59.00]
    bn254_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [59.00, 154.03]
    bn254_G2::fixed_base_exp_window_table.push_back(59);
    // window 6 is unbeaten in [154.03, 334.25]
    bn254_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [334.25, 742.58]
    bn254_G2::fixed_base_exp_window_table.push_back(334);
    // window 8 is unbeaten in [742.58, 2034.40]
    bn254_G2::fixed_base_exp_window_table.push_back(743);
    // window 9 is unbeaten in [2034.40, 4987.56]
    bn254_G2::fixed_base_exp_window_table.push_back(2034);
    // window 10 is unbeaten in [4987.56, 8888.27]
    bn254_G2::fixed_base_exp_window_table.push_back(4988);
    // window 11 is unbeaten in [8888.27, 26271.13]
    bn254_G2::fixed_base_exp_window_table.push_back(8888);
    // window 12 is unbeaten in [26271.13, 39768.20]
    bn254_G2::fixed_base_exp_window_table.push_back(26271);
    // window 13 is unbeaten in [39768.20, 106275.75]
    bn254_G2::fixed_base_exp_window_table.push_back(39768);
    // window 14 is unbeaten in [106275.75, 141703.40]
    bn254_G2::fixed_base_exp_window_table.push_back(106276);
    // window 15 is unbeaten in [141703.40, 462422.97]
    bn254_G2::fixed_base_exp_window_table.push_back(141703);
    // window 16 is unbeaten in [462422.97, 926871.84]
    bn254_G2::fixed_base_exp_window_table.push_back(462423);
    // window 17 is unbeaten in [926871.84, 4873049.17]
    bn254_G2::fixed_base_exp_window_table.push_back(926872);
    // window 18 is never the best
    bn254_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4873049.17, 5706707.88]
    bn254_G2::fixed_base_exp_window_table.push_back(4873049);
    // window 20 is unbeaten in [5706707.88, 31673814.95]
    bn254_G2::fixed_base_exp_window_table.push_back(5706708);
    // window 21 is never the best
    bn254_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [31673814.95, inf]
    bn254_G2::fixed_base_exp_window_table.push_back(31673815);

    /* choice of pairing */

	bn254_ate_loop_count = bigint<bn254_q_limbs>("29793968203157093288");
	bn254_ate_is_loop_count_neg = false;
	bn254_final_exponent_z = bigint<bn254_q_limbs>("4965661367192848881");
	bn254_final_exponent_is_z_neg = false;

}
} // libff
