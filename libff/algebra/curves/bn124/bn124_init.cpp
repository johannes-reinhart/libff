#include <libff/algebra/curves/bn124/bn124_g1.hpp>
#include <libff/algebra/curves/bn124/bn124_g2.hpp>
#include <libff/algebra/curves/bn124/bn124_init.hpp>

namespace libff {

bn124_Fq bn124_coeff_b;
bn124_Fq2 bn124_twist;
bn124_Fq2 bn124_twist_coeff_b;
bn124_Fq bn124_twist_mul_by_b_c0;
bn124_Fq bn124_twist_mul_by_b_c1;
bn124_Fq2 bn124_twist_mul_by_q_X;
bn124_Fq2 bn124_twist_mul_by_q_Y;

bigint<bn124_q_limbs> bn124_ate_loop_count;
bool bn124_ate_is_loop_count_neg;
bigint<bn124_q_limbs> bn124_final_exponent_z;
bool bn124_final_exponent_is_z_neg;

void init_bn124_params()
{
	init_bn124_fields();

    /* choice of short Weierstrass curve and its twist */

	bn124_coeff_b = bn124_Fq("3");
	bn124_twist = bn124_Fq2(bn124_Fq("5"), bn124_Fq("1"));
	bn124_twist_coeff_b = bn124_coeff_b * bn124_twist.inverse();
	bn124_twist_mul_by_b_c0 = bn124_coeff_b * bn124_Fq2::non_residue;
	bn124_twist_mul_by_b_c1 = bn124_coeff_b * bn124_Fq2::non_residue;
	bn124_twist_mul_by_q_X = bn124_Fq2(bn124_Fq("11057421525851115170548223099100677627"),bn124_Fq("5501413929240230020448532534432497469"));
	bn124_twist_mul_by_q_Y = bn124_Fq2(bn124_Fq("11217235539836039461521625978614135522"),bn124_Fq("5085777724801701118551070705763167893"));
	
    /* choice of group G1 */

	bn124_G1::G1_zero = bn124_G1(bn124_Fq::zero(),   
                                                   bn124_Fq::one(),
                                                   bn124_Fq::zero());
	bn124_G1::G1_one = bn124_G1(bn124_Fq("1"),
                                                  bn124_Fq("2"),
                                                  bn124_Fq::one());
    bn124_G1::initialized = true;

    // Cofactor
    bn124_G1::h = bigint<bn124_G1::h_limbs>("1");

    // TODO
	bn124_G1::wnaf_window_table.resize(0);
    bn124_G1::wnaf_window_table.push_back(11);
    bn124_G1::wnaf_window_table.push_back(24);
    bn124_G1::wnaf_window_table.push_back(60);
    bn124_G1::wnaf_window_table.push_back(127);

    bn124_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.99]
    bn124_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.99, 10.99]
    bn124_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.99, 32.29]
    bn124_G1::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [32.29, 55.23]
    bn124_G1::fixed_base_exp_window_table.push_back(32);
    // window 5 is unbeaten in [55.23, 162.03]
    bn124_G1::fixed_base_exp_window_table.push_back(55);
    // window 6 is unbeaten in [162.03, 360.15]
    bn124_G1::fixed_base_exp_window_table.push_back(162);
    // window 7 is unbeaten in [360.15, 815.44]
    bn124_G1::fixed_base_exp_window_table.push_back(360);
    // window 8 is unbeaten in [815.44, 2373.07]
    bn124_G1::fixed_base_exp_window_table.push_back(815);
    // window 9 is unbeaten in [2373.07, 6977.75]
    bn124_G1::fixed_base_exp_window_table.push_back(2373);
    // window 10 is unbeaten in [6977.75, 7122.23]
    bn124_G1::fixed_base_exp_window_table.push_back(6978);
    // window 11 is unbeaten in [7122.23, 57818.46]
    bn124_G1::fixed_base_exp_window_table.push_back(7122);
    // window 12 is never the best
    bn124_G1::fixed_base_exp_window_table.push_back(0);
    // window 13 is unbeaten in [57818.46, 169679.14]
    bn124_G1::fixed_base_exp_window_table.push_back(57818);
    // window 14 is never the best
    bn124_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [169679.14, 439758.91]
    bn124_G1::fixed_base_exp_window_table.push_back(169679);
    // window 16 is unbeaten in [439758.91, 936073.41]
    bn124_G1::fixed_base_exp_window_table.push_back(439759);
    // window 17 is unbeaten in [936073.41, 4666554.74]
    bn124_G1::fixed_base_exp_window_table.push_back(936073);
    // window 18 is never the best
    bn124_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4666554.74, 7580404.42]
    bn124_G1::fixed_base_exp_window_table.push_back(4666555);
    // window 20 is unbeaten in [7580404.42, 34552892.20]
    bn124_G1::fixed_base_exp_window_table.push_back(7580404);
    // window 21 is never the best
    bn124_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [34552892.20, inf]
    bn124_G1::fixed_base_exp_window_table.push_back(34552892);

    /* choice of group G2 */

	bn124_G2::G2_zero = bn124_G2(bn124_Fq2::zero(),
                                 bn124_Fq2::one(),
                                 bn124_Fq2::zero());
	bn124_G2::G2_one = bn124_G2(bn124_Fq2(bn124_Fq("1482181262851106128340730184604596283"),
                                                                     bn124_Fq("10066821713171670698249225372320188548")),
                                                  bn124_Fq2(bn124_Fq("1586481144356933234099397925323010562"),
                                                                     bn124_Fq("387806682866053068498307761784234099")),
                                                  bn124_Fq2::one());

    // Cofactor
    bn124_G2::h = bigint<bn124_G2::h_limbs>("17000133324792832067142141520207542925");

    // TODO
	bn124_G2::wnaf_window_table.resize(0);
    bn124_G2::wnaf_window_table.push_back(5);
    bn124_G2::wnaf_window_table.push_back(15);
    bn124_G2::wnaf_window_table.push_back(39);
    bn124_G2::wnaf_window_table.push_back(109);

    bn124_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.10]
    bn124_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.10, 10.43]
    bn124_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.43, 25.28]
    bn124_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.28, 59.00]
    bn124_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [59.00, 154.03]
    bn124_G2::fixed_base_exp_window_table.push_back(59);
    // window 6 is unbeaten in [154.03, 334.25]
    bn124_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [334.25, 742.58]
    bn124_G2::fixed_base_exp_window_table.push_back(334);
    // window 8 is unbeaten in [742.58, 2034.40]
    bn124_G2::fixed_base_exp_window_table.push_back(743);
    // window 9 is unbeaten in [2034.40, 4987.56]
    bn124_G2::fixed_base_exp_window_table.push_back(2034);
    // window 10 is unbeaten in [4987.56, 8888.27]
    bn124_G2::fixed_base_exp_window_table.push_back(4988);
    // window 11 is unbeaten in [8888.27, 26271.13]
    bn124_G2::fixed_base_exp_window_table.push_back(8888);
    // window 12 is unbeaten in [26271.13, 39768.20]
    bn124_G2::fixed_base_exp_window_table.push_back(26271);
    // window 13 is unbeaten in [39768.20, 106275.75]
    bn124_G2::fixed_base_exp_window_table.push_back(39768);
    // window 14 is unbeaten in [106275.75, 141703.40]
    bn124_G2::fixed_base_exp_window_table.push_back(106276);
    // window 15 is unbeaten in [141703.40, 462422.97]
    bn124_G2::fixed_base_exp_window_table.push_back(141703);
    // window 16 is unbeaten in [462422.97, 926871.84]
    bn124_G2::fixed_base_exp_window_table.push_back(462423);
    // window 17 is unbeaten in [926871.84, 4873049.17]
    bn124_G2::fixed_base_exp_window_table.push_back(926872);
    // window 18 is never the best
    bn124_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4873049.17, 5706707.88]
    bn124_G2::fixed_base_exp_window_table.push_back(4873049);
    // window 20 is unbeaten in [5706707.88, 31673814.95]
    bn124_G2::fixed_base_exp_window_table.push_back(5706708);
    // window 21 is never the best
    bn124_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [31673814.95, inf]
    bn124_G2::fixed_base_exp_window_table.push_back(31673815);

    /* choice of pairing */

	bn124_ate_loop_count = bigint<bn124_q_limbs>("4973804456");
	bn124_ate_is_loop_count_neg = false;
	bn124_final_exponent_z = bigint<bn124_q_limbs>("828967409");
	bn124_final_exponent_is_z_neg = false;

}
} // libff
