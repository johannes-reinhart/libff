#include <libff/algebra/curves/edwards58/edwards58_g1.hpp>
#include <libff/algebra/curves/edwards58/edwards58_g2.hpp>
#include <libff/algebra/curves/edwards58/edwards58_init.hpp>

namespace libff {

edwards58_Fq3 edwards58_twist;
edwards58_Fq3 edwards58_twist_coeff_a;
edwards58_Fq3 edwards58_twist_coeff_d;
edwards58_Fq edwards58_twist_mul_by_a_c0;
edwards58_Fq edwards58_twist_mul_by_a_c1;
edwards58_Fq edwards58_twist_mul_by_a_c2;
edwards58_Fq edwards58_twist_mul_by_d_c0;
edwards58_Fq edwards58_twist_mul_by_d_c1;
edwards58_Fq edwards58_twist_mul_by_d_c2;
edwards58_Fq edwards58_twist_mul_by_q_Y;
edwards58_Fq edwards58_twist_mul_by_q_Z;

bigint<edwards58_q_limbs> edwards58_ate_loop_count;
bigint<edwards58_q_limbs> edwards58_final_exponent_last_chunk_abs_of_w0;
bool edwards58_final_exponent_last_chunk_is_w0_neg;
bigint<edwards58_q_limbs> edwards58_final_exponent_last_chunk_w1;

void init_edwards58_params()
{
    init_edwards58_fields();

    /* choice of Edwards curve and its twist */

    edwards58_G1::coeff_a = edwards58_Fq("5");
    edwards58_G1::coeff_d = edwards58_Fq("579073710274753001");
    edwards58_twist = edwards58_Fq3(edwards58_Fq::zero(), edwards58_Fq::one(), edwards58_Fq::zero());
    edwards58_twist_coeff_a = edwards58_G1::coeff_a * edwards58_twist;
    edwards58_twist_coeff_d = edwards58_G1::coeff_d * edwards58_twist;
    edwards58_twist_mul_by_a_c0 = edwards58_G1::coeff_a * edwards58_Fq3::non_residue;
    edwards58_twist_mul_by_a_c1 = edwards58_G1::coeff_a;
    edwards58_twist_mul_by_a_c2 = edwards58_G1::coeff_a;
    edwards58_twist_mul_by_d_c0 = edwards58_G1::coeff_d * edwards58_Fq3::non_residue;
    edwards58_twist_mul_by_d_c1 = edwards58_G1::coeff_d;
    edwards58_twist_mul_by_d_c2 = edwards58_G1::coeff_d;
    edwards58_twist_mul_by_q_Y = edwards58_Fq("126784148197437871");
    edwards58_twist_mul_by_q_Z = edwards58_Fq("126784148197437871");

    /* choice of group G1 */
    edwards58_G1::G1_zero = edwards58_G1(edwards58_Fq::zero(),
                                     edwards58_Fq::one());
    edwards58_G1::G1_one = edwards58_G1(edwards58_Fq("135119008168998470"),
                                    edwards58_Fq("793956801732934748"));
    edwards58_G1::initialized = true;

    // TODO
    edwards58_G1::wnaf_window_table.resize(0);
    edwards58_G1::wnaf_window_table.push_back(9);
    edwards58_G1::wnaf_window_table.push_back(14);
    edwards58_G1::wnaf_window_table.push_back(24);
    edwards58_G1::wnaf_window_table.push_back(117);

    edwards58_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.10]
    edwards58_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.10, 9.69]
    edwards58_G1::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [9.69, 25.21]
    edwards58_G1::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.21, 60.00]
    edwards58_G1::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.00, 149.33]
    edwards58_G1::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [149.33, 369.61]
    edwards58_G1::fixed_base_exp_window_table.push_back(149);
    // window 7 is unbeaten in [369.61, 849.07]
    edwards58_G1::fixed_base_exp_window_table.push_back(370);
    // window 8 is unbeaten in [849.07, 1764.94]
    edwards58_G1::fixed_base_exp_window_table.push_back(849);
    // window 9 is unbeaten in [1764.94, 4429.59]
    edwards58_G1::fixed_base_exp_window_table.push_back(1765);
    // window 10 is unbeaten in [4429.59, 13388.78]
    edwards58_G1::fixed_base_exp_window_table.push_back(4430);
    // window 11 is unbeaten in [13388.78, 15368.00]
    edwards58_G1::fixed_base_exp_window_table.push_back(13389);
    // window 12 is unbeaten in [15368.00, 74912.07]
    edwards58_G1::fixed_base_exp_window_table.push_back(15368);
    // window 13 is unbeaten in [74912.07, 438107.20]
    edwards58_G1::fixed_base_exp_window_table.push_back(74912);
    // window 14 is never the best
    edwards58_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [438107.20, 1045626.18]
    edwards58_G1::fixed_base_exp_window_table.push_back(438107);
    // window 16 is never the best
    edwards58_G1::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [1045626.18, 1577434.48]
    edwards58_G1::fixed_base_exp_window_table.push_back(1045626);
    // window 18 is unbeaten in [1577434.48, 17350594.23]
    edwards58_G1::fixed_base_exp_window_table.push_back(1577434);
    // window 19 is never the best
    edwards58_G1::fixed_base_exp_window_table.push_back(0);
    // window 20 is never the best
    edwards58_G1::fixed_base_exp_window_table.push_back(0);
    // window 21 is unbeaten in [17350594.23, inf]
    edwards58_G1::fixed_base_exp_window_table.push_back(17350594);
    // window 22 is never the best
    edwards58_G1::fixed_base_exp_window_table.push_back(0);

    /* choice of group G2 */

    edwards58_G2::G2_zero = edwards58_G2(edwards58_Fq3::zero(),
                                     edwards58_Fq3::one());
    edwards58_G2::G2_one = edwards58_G2(edwards58_Fq3(edwards58_Fq("686452412557329133"),
                                                edwards58_Fq("485520559859571402"),
                                                edwards58_Fq("238985303084351179")),
                                    edwards58_Fq3(edwards58_Fq("114409345708660176"),
                                                edwards58_Fq("241586179836920758"),
                                                edwards58_Fq("154947719077520267")));
    edwards58_G2::initialized = true;

    // TODO
    edwards58_G2::wnaf_window_table.resize(0);
    edwards58_G2::wnaf_window_table.push_back(6);
    edwards58_G2::wnaf_window_table.push_back(12);
    edwards58_G2::wnaf_window_table.push_back(42);
    edwards58_G2::wnaf_window_table.push_back(97);

    edwards58_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.74]
    edwards58_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.74, 10.67]
    edwards58_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.67, 25.53]
    edwards58_G2::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [25.53, 60.67]
    edwards58_G2::fixed_base_exp_window_table.push_back(26);
    // window 5 is unbeaten in [60.67, 145.77]
    edwards58_G2::fixed_base_exp_window_table.push_back(61);
    // window 6 is unbeaten in [145.77, 356.76]
    edwards58_G2::fixed_base_exp_window_table.push_back(146);
    // window 7 is unbeaten in [356.76, 823.08]
    edwards58_G2::fixed_base_exp_window_table.push_back(357);
    // window 8 is unbeaten in [823.08, 1589.45]
    edwards58_G2::fixed_base_exp_window_table.push_back(823);
    // window 9 is unbeaten in [1589.45, 4135.70]
    edwards58_G2::fixed_base_exp_window_table.push_back(1589);
    // window 10 is unbeaten in [4135.70, 14297.74]
    edwards58_G2::fixed_base_exp_window_table.push_back(4136);
    // window 11 is unbeaten in [14297.74, 16744.85]
    edwards58_G2::fixed_base_exp_window_table.push_back(14298);
    // window 12 is unbeaten in [16744.85, 51768.98]
    edwards58_G2::fixed_base_exp_window_table.push_back(16745);
    // window 13 is unbeaten in [51768.98, 99811.01]
    edwards58_G2::fixed_base_exp_window_table.push_back(51769);
    // window 14 is unbeaten in [99811.01, 193306.72]
    edwards58_G2::fixed_base_exp_window_table.push_back(99811);
    // window 15 is unbeaten in [193306.72, 907184.68]
    edwards58_G2::fixed_base_exp_window_table.push_back(193307);
    // window 16 is never the best
    edwards58_G2::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [907184.68, 1389682.59]
    edwards58_G2::fixed_base_exp_window_table.push_back(907185);
    // window 18 is unbeaten in [1389682.59, 6752695.74]
    edwards58_G2::fixed_base_exp_window_table.push_back(1389683);
    // window 19 is never the best
    edwards58_G2::fixed_base_exp_window_table.push_back(0);
    // window 20 is unbeaten in [6752695.74, 193642894.51]
    edwards58_G2::fixed_base_exp_window_table.push_back(6752696);
    // window 21 is unbeaten in [193642894.51, 226760202.29]
    edwards58_G2::fixed_base_exp_window_table.push_back(193642895);
    // window 22 is unbeaten in [226760202.29, inf]
    edwards58_G2::fixed_base_exp_window_table.push_back(226760202);

    /* pairing parameters */

    edwards58_ate_loop_count = bigint<edwards58_q_limbs>("1656225795");
    edwards58_final_exponent_last_chunk_abs_of_w0 = bigint<edwards58_q_limbs>("6624903171");
    edwards58_final_exponent_last_chunk_is_w0_neg = true;
    edwards58_final_exponent_last_chunk_w1 = bigint<edwards58_q_limbs>("4");


}
} // libff
