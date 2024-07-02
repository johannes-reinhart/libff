#include <libff/algebra/curves/bn183/bn183_g1.hpp>
#include <libff/algebra/curves/bn183/bn183_g2.hpp>
#include <libff/algebra/curves/bn183/bn183_init.hpp>

namespace libff {

bn183_Fq bn183_coeff_b;
bn183_Fq2 bn183_twist;
bn183_Fq2 bn183_twist_coeff_b;
bn183_Fq bn183_twist_mul_by_b_c0;
bn183_Fq bn183_twist_mul_by_b_c1;
bn183_Fq2 bn183_twist_mul_by_q_X;
bn183_Fq2 bn183_twist_mul_by_q_Y;

bigint<bn183_q_limbs> bn183_ate_loop_count;
bool bn183_ate_is_loop_count_neg;
bigint<bn183_q_limbs> bn183_final_exponent_z;
bool bn183_final_exponent_is_z_neg;

void init_bn183_params()
{
	init_bn183_fields();

    /* choice of short Weierstrass curve and its twist */

	bn183_coeff_b = bn183_Fq("3");
	bn183_twist = bn183_Fq2(bn183_Fq("2"), bn183_Fq("1"));
	bn183_twist_coeff_b = bn183_coeff_b * bn183_twist.inverse();
	bn183_twist_mul_by_b_c0 = bn183_coeff_b * bn183_Fq2::non_residue;
	bn183_twist_mul_by_b_c1 = bn183_coeff_b * bn183_Fq2::non_residue;
	bn183_twist_mul_by_q_X = bn183_Fq2(bn183_Fq("2514746150702782023387381073984499065469793711673933027"),bn183_Fq("3170995372818140140778075409082321578161037698420449167"));
	bn183_twist_mul_by_q_Y = bn183_Fq2(bn183_Fq("2785947042921859904240819706842983078799881943539725802"),bn183_Fq("5571894085843719808481639413685966157599763887079451604"));
	
    /* choice of group G1 */

	bn183_G1::G1_zero = bn183_G1(bn183_Fq::zero(),   
                                                   bn183_Fq::one(),
                                                   bn183_Fq::zero());
	bn183_G1::G1_one = bn183_G1(bn183_Fq("1"),
                                                  bn183_Fq("2"),
                                                  bn183_Fq::one());
    bn183_G1::initialized = true;

    // Cofactor
    bn183_G1::h = bigint<bn183_G1::h_limbs>("1");

    // TODO
	bn183_G1::wnaf_window_table.resize(0);
    bn183_G1::wnaf_window_table.push_back(11);
    bn183_G1::wnaf_window_table.push_back(24);
    bn183_G1::wnaf_window_table.push_back(60);
    bn183_G1::wnaf_window_table.push_back(127);

    bn183_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.99]
    bn183_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.99, 10.99]
    bn183_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.99, 32.29]
    bn183_G1::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [32.29, 55.23]
    bn183_G1::fixed_base_exp_window_table.push_back(32);
    // window 5 is unbeaten in [55.23, 162.03]
    bn183_G1::fixed_base_exp_window_table.push_back(55);
    // window 6 is unbeaten in [162.03, 360.15]
    bn183_G1::fixed_base_exp_window_table.push_back(162);
    // window 7 is unbeaten in [360.15, 815.44]
    bn183_G1::fixed_base_exp_window_table.push_back(360);
    // window 8 is unbeaten in [815.44, 2373.07]
    bn183_G1::fixed_base_exp_window_table.push_back(815);
    // window 9 is unbeaten in [2373.07, 6977.75]
    bn183_G1::fixed_base_exp_window_table.push_back(2373);
    // window 10 is unbeaten in [6977.75, 7122.23]
    bn183_G1::fixed_base_exp_window_table.push_back(6978);
    // window 11 is unbeaten in [7122.23, 57818.46]
    bn183_G1::fixed_base_exp_window_table.push_back(7122);
    // window 12 is never the best
    bn183_G1::fixed_base_exp_window_table.push_back(0);
    // window 13 is unbeaten in [57818.46, 169679.14]
    bn183_G1::fixed_base_exp_window_table.push_back(57818);
    // window 14 is never the best
    bn183_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [169679.14, 439758.91]
    bn183_G1::fixed_base_exp_window_table.push_back(169679);
    // window 16 is unbeaten in [439758.91, 936073.41]
    bn183_G1::fixed_base_exp_window_table.push_back(439759);
    // window 17 is unbeaten in [936073.41, 4666554.74]
    bn183_G1::fixed_base_exp_window_table.push_back(936073);
    // window 18 is never the best
    bn183_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4666554.74, 7580404.42]
    bn183_G1::fixed_base_exp_window_table.push_back(4666555);
    // window 20 is unbeaten in [7580404.42, 34552892.20]
    bn183_G1::fixed_base_exp_window_table.push_back(7580404);
    // window 21 is never the best
    bn183_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [34552892.20, inf]
    bn183_G1::fixed_base_exp_window_table.push_back(34552892);

    /* choice of group G2 */

	bn183_G2::G2_zero = bn183_G2(bn183_Fq2::zero(),
                                 bn183_Fq2::one(),
                                 bn183_Fq2::zero());
	bn183_G2::G2_one = bn183_G2(bn183_Fq2(bn183_Fq("721678245984884345543303064655407316345703950658223194"),
                                                                     bn183_Fq("2575758008344461810840893514358025270228962025833591182")),
                                                  bn183_Fq2(bn183_Fq("1908382470610824589133495396746666294729025620843352632"),
                                                                     bn183_Fq("3544812708057231864739467383540395526301995413400188899")),
                                                  bn183_Fq2::one());

    // Cofactor
    bn183_G2::h = bigint<bn183_G2::h_limbs>("6804759748846355405830582791228219856011899129244023437");

    // TODO
	bn183_G2::wnaf_window_table.resize(0);
    bn183_G2::wnaf_window_table.push_back(5);
    bn183_G2::wnaf_window_table.push_back(15);
    bn183_G2::wnaf_window_table.push_back(39);
    bn183_G2::wnaf_window_table.push_back(109);

    bn183_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.10]
    bn183_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.10, 10.43]
    bn183_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.43, 25.28]
    bn183_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.28, 59.00]
    bn183_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [59.00, 154.03]
    bn183_G2::fixed_base_exp_window_table.push_back(59);
    // window 6 is unbeaten in [154.03, 334.25]
    bn183_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [334.25, 742.58]
    bn183_G2::fixed_base_exp_window_table.push_back(334);
    // window 8 is unbeaten in [742.58, 2034.40]
    bn183_G2::fixed_base_exp_window_table.push_back(743);
    // window 9 is unbeaten in [2034.40, 4987.56]
    bn183_G2::fixed_base_exp_window_table.push_back(2034);
    // window 10 is unbeaten in [4987.56, 8888.27]
    bn183_G2::fixed_base_exp_window_table.push_back(4988);
    // window 11 is unbeaten in [8888.27, 26271.13]
    bn183_G2::fixed_base_exp_window_table.push_back(8888);
    // window 12 is unbeaten in [26271.13, 39768.20]
    bn183_G2::fixed_base_exp_window_table.push_back(26271);
    // window 13 is unbeaten in [39768.20, 106275.75]
    bn183_G2::fixed_base_exp_window_table.push_back(39768);
    // window 14 is unbeaten in [106275.75, 141703.40]
    bn183_G2::fixed_base_exp_window_table.push_back(106276);
    // window 15 is unbeaten in [141703.40, 462422.97]
    bn183_G2::fixed_base_exp_window_table.push_back(141703);
    // window 16 is unbeaten in [462422.97, 926871.84]
    bn183_G2::fixed_base_exp_window_table.push_back(462423);
    // window 17 is unbeaten in [926871.84, 4873049.17]
    bn183_G2::fixed_base_exp_window_table.push_back(926872);
    // window 18 is never the best
    bn183_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4873049.17, 5706707.88]
    bn183_G2::fixed_base_exp_window_table.push_back(4873049);
    // window 20 is unbeaten in [5706707.88, 31673814.95]
    bn183_G2::fixed_base_exp_window_table.push_back(5706708);
    // window 21 is never the best
    bn183_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [31673814.95, inf]
    bn183_G2::fixed_base_exp_window_table.push_back(31673815);

    /* choice of pairing */

	bn183_ate_loop_count = bigint<bn183_q_limbs>("125106197511080");
	bn183_ate_is_loop_count_neg = false;
	bn183_final_exponent_z = bigint<bn183_q_limbs>("20851032918513");
	bn183_final_exponent_is_z_neg = false;

}
} // libff
