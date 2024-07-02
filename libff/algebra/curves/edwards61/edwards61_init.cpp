#include <libff/algebra/curves/edwards61/edwards61_g1.hpp>
#include <libff/algebra/curves/edwards61/edwards61_g2.hpp>
#include <libff/algebra/curves/edwards61/edwards61_init.hpp>

namespace libff {

edwards61_Fq3 edwards61_twist;
edwards61_Fq3 edwards61_twist_coeff_a;
edwards61_Fq3 edwards61_twist_coeff_d;
edwards61_Fq edwards61_twist_mul_by_a_c0;
edwards61_Fq edwards61_twist_mul_by_a_c1;
edwards61_Fq edwards61_twist_mul_by_a_c2;
edwards61_Fq edwards61_twist_mul_by_d_c0;
edwards61_Fq edwards61_twist_mul_by_d_c1;
edwards61_Fq edwards61_twist_mul_by_d_c2;
edwards61_Fq edwards61_twist_mul_by_q_Y;
edwards61_Fq edwards61_twist_mul_by_q_Z;

bigint<edwards61_q_limbs> edwards61_ate_loop_count;
bigint<edwards61_q_limbs> edwards61_final_exponent_last_chunk_abs_of_w0;
bool edwards61_final_exponent_last_chunk_is_w0_neg;
bigint<edwards61_q_limbs> edwards61_final_exponent_last_chunk_w1;

void init_edwards61_params()
{
    init_edwards61_fields();

    /* choice of Edwards curve and its twist */

    edwards61_G1::coeff_a = edwards61_Fq("5");
    edwards61_G1::coeff_d = edwards61_Fq("2154376507894293213");
    edwards61_twist = edwards61_Fq3(edwards61_Fq::zero(), edwards61_Fq::one(), edwards61_Fq::zero());
    edwards61_twist_coeff_a = edwards61_G1::coeff_a * edwards61_twist;
    edwards61_twist_coeff_d = edwards61_G1::coeff_d * edwards61_twist;
    edwards61_twist_mul_by_a_c0 = edwards61_G1::coeff_a * edwards61_Fq3::non_residue;
    edwards61_twist_mul_by_a_c1 = edwards61_G1::coeff_a;
    edwards61_twist_mul_by_a_c2 = edwards61_G1::coeff_a;
    edwards61_twist_mul_by_d_c0 = edwards61_G1::coeff_d * edwards61_Fq3::non_residue;
    edwards61_twist_mul_by_d_c1 = edwards61_G1::coeff_d;
    edwards61_twist_mul_by_d_c2 = edwards61_G1::coeff_d;
    edwards61_twist_mul_by_q_Y = edwards61_Fq("875438387576607506");
    edwards61_twist_mul_by_q_Z = edwards61_Fq("875438387576607506");

    /* choice of group G1 */
    edwards61_G1::G1_zero = edwards61_G1(edwards61_Fq::zero(),
                                     edwards61_Fq::one());
    edwards61_G1::G1_one = edwards61_G1(edwards61_Fq("2564941438079622667"),
                                    edwards61_Fq("268572741021024558"));
    edwards61_G1::initialized = true;

    // TODO
    edwards61_G1::wnaf_window_table.resize(0);
    edwards61_G1::wnaf_window_table.push_back(9);
    edwards61_G1::wnaf_window_table.push_back(14);
    edwards61_G1::wnaf_window_table.push_back(24);
    edwards61_G1::wnaf_window_table.push_back(117);

    edwards61_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.10]
    edwards61_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.10, 9.69]
    edwards61_G1::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [9.69, 25.21]
    edwards61_G1::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.21, 60.00]
    edwards61_G1::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.00, 149.33]
    edwards61_G1::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [149.33, 369.61]
    edwards61_G1::fixed_base_exp_window_table.push_back(149);
    // window 7 is unbeaten in [369.61, 849.07]
    edwards61_G1::fixed_base_exp_window_table.push_back(370);
    // window 8 is unbeaten in [849.07, 1764.94]
    edwards61_G1::fixed_base_exp_window_table.push_back(849);
    // window 9 is unbeaten in [1764.94, 4429.59]
    edwards61_G1::fixed_base_exp_window_table.push_back(1765);
    // window 10 is unbeaten in [4429.59, 13388.78]
    edwards61_G1::fixed_base_exp_window_table.push_back(4430);
    // window 11 is unbeaten in [13388.78, 15368.00]
    edwards61_G1::fixed_base_exp_window_table.push_back(13389);
    // window 12 is unbeaten in [15368.00, 74912.07]
    edwards61_G1::fixed_base_exp_window_table.push_back(15368);
    // window 13 is unbeaten in [74912.07, 438107.20]
    edwards61_G1::fixed_base_exp_window_table.push_back(74912);
    // window 14 is never the best
    edwards61_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [438107.20, 1045626.18]
    edwards61_G1::fixed_base_exp_window_table.push_back(438107);
    // window 16 is never the best
    edwards61_G1::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [1045626.18, 1577434.48]
    edwards61_G1::fixed_base_exp_window_table.push_back(1045626);
    // window 18 is unbeaten in [1577434.48, 17350594.23]
    edwards61_G1::fixed_base_exp_window_table.push_back(1577434);
    // window 19 is never the best
    edwards61_G1::fixed_base_exp_window_table.push_back(0);
    // window 20 is never the best
    edwards61_G1::fixed_base_exp_window_table.push_back(0);
    // window 21 is unbeaten in [17350594.23, inf]
    edwards61_G1::fixed_base_exp_window_table.push_back(17350594);
    // window 22 is never the best
    edwards61_G1::fixed_base_exp_window_table.push_back(0);

    /* choice of group G2 */

    edwards61_G2::G2_zero = edwards61_G2(edwards61_Fq3::zero(),
                                     edwards61_Fq3::one());
    edwards61_G2::G2_one = edwards61_G2(edwards61_Fq3(edwards61_Fq("2867310054428652554"),
                                                edwards61_Fq("3210449111018192129"),
                                                edwards61_Fq("2899260707215196651")),
                                    edwards61_Fq3(edwards61_Fq("2683722972417512066"),
                                                edwards61_Fq("4524248731112963787"),
                                                edwards61_Fq("1809298747241172860")));
    edwards61_G2::initialized = true;

    // TODO
    edwards61_G2::wnaf_window_table.resize(0);
    edwards61_G2::wnaf_window_table.push_back(6);
    edwards61_G2::wnaf_window_table.push_back(12);
    edwards61_G2::wnaf_window_table.push_back(42);
    edwards61_G2::wnaf_window_table.push_back(97);

    edwards61_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.74]
    edwards61_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.74, 10.67]
    edwards61_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.67, 25.53]
    edwards61_G2::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [25.53, 60.67]
    edwards61_G2::fixed_base_exp_window_table.push_back(26);
    // window 5 is unbeaten in [60.67, 145.77]
    edwards61_G2::fixed_base_exp_window_table.push_back(61);
    // window 6 is unbeaten in [145.77, 356.76]
    edwards61_G2::fixed_base_exp_window_table.push_back(146);
    // window 7 is unbeaten in [356.76, 823.08]
    edwards61_G2::fixed_base_exp_window_table.push_back(357);
    // window 8 is unbeaten in [823.08, 1589.45]
    edwards61_G2::fixed_base_exp_window_table.push_back(823);
    // window 9 is unbeaten in [1589.45, 4135.70]
    edwards61_G2::fixed_base_exp_window_table.push_back(1589);
    // window 10 is unbeaten in [4135.70, 14297.74]
    edwards61_G2::fixed_base_exp_window_table.push_back(4136);
    // window 11 is unbeaten in [14297.74, 16744.85]
    edwards61_G2::fixed_base_exp_window_table.push_back(14298);
    // window 12 is unbeaten in [16744.85, 51768.98]
    edwards61_G2::fixed_base_exp_window_table.push_back(16745);
    // window 13 is unbeaten in [51768.98, 99811.01]
    edwards61_G2::fixed_base_exp_window_table.push_back(51769);
    // window 14 is unbeaten in [99811.01, 193306.72]
    edwards61_G2::fixed_base_exp_window_table.push_back(99811);
    // window 15 is unbeaten in [193306.72, 907184.68]
    edwards61_G2::fixed_base_exp_window_table.push_back(193307);
    // window 16 is never the best
    edwards61_G2::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [907184.68, 1389682.59]
    edwards61_G2::fixed_base_exp_window_table.push_back(907185);
    // window 18 is unbeaten in [1389682.59, 6752695.74]
    edwards61_G2::fixed_base_exp_window_table.push_back(1389683);
    // window 19 is never the best
    edwards61_G2::fixed_base_exp_window_table.push_back(0);
    // window 20 is unbeaten in [6752695.74, 193642894.51]
    edwards61_G2::fixed_base_exp_window_table.push_back(6752696);
    // window 21 is unbeaten in [193642894.51, 226760202.29]
    edwards61_G2::fixed_base_exp_window_table.push_back(193642895);
    // window 22 is unbeaten in [226760202.29, inf]
    edwards61_G2::fixed_base_exp_window_table.push_back(226760202);

    /* pairing parameters */

    edwards61_ate_loop_count = bigint<edwards61_q_limbs>("4007657475");
    edwards61_final_exponent_last_chunk_abs_of_w0 = bigint<edwards61_q_limbs>("16030629891");
    edwards61_final_exponent_last_chunk_is_w0_neg = true;
    edwards61_final_exponent_last_chunk_w1 = bigint<edwards61_q_limbs>("4");


}
} // libff
