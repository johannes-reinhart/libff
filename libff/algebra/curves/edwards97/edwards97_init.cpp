#include <libff/algebra/curves/edwards97/edwards97_g1.hpp>
#include <libff/algebra/curves/edwards97/edwards97_g2.hpp>
#include <libff/algebra/curves/edwards97/edwards97_init.hpp>

namespace libff {

edwards97_Fq3 edwards97_twist;
edwards97_Fq3 edwards97_twist_coeff_a;
edwards97_Fq3 edwards97_twist_coeff_d;
edwards97_Fq edwards97_twist_mul_by_a_c0;
edwards97_Fq edwards97_twist_mul_by_a_c1;
edwards97_Fq edwards97_twist_mul_by_a_c2;
edwards97_Fq edwards97_twist_mul_by_d_c0;
edwards97_Fq edwards97_twist_mul_by_d_c1;
edwards97_Fq edwards97_twist_mul_by_d_c2;
edwards97_Fq edwards97_twist_mul_by_q_Y;
edwards97_Fq edwards97_twist_mul_by_q_Z;

bigint<edwards97_q_limbs> edwards97_ate_loop_count;
bigint<edwards97_q_limbs> edwards97_final_exponent_last_chunk_abs_of_w0;
bool edwards97_final_exponent_last_chunk_is_w0_neg;
bigint<edwards97_q_limbs> edwards97_final_exponent_last_chunk_w1;

void init_edwards97_params()
{
    init_edwards97_fields();

    /* choice of Edwards curve and its twist */

    edwards97_G1::coeff_a = edwards97_Fq("5");
    edwards97_G1::coeff_d = edwards97_Fq("482996825047815773983380486779");
    edwards97_twist = edwards97_Fq3(edwards97_Fq::zero(), edwards97_Fq::one(), edwards97_Fq::zero());
    edwards97_twist_coeff_a = edwards97_G1::coeff_a * edwards97_twist;
    edwards97_twist_coeff_d = edwards97_G1::coeff_d * edwards97_twist;
    edwards97_twist_mul_by_a_c0 = edwards97_G1::coeff_a * edwards97_Fq3::non_residue;
    edwards97_twist_mul_by_a_c1 = edwards97_G1::coeff_a;
    edwards97_twist_mul_by_a_c2 = edwards97_G1::coeff_a;
    edwards97_twist_mul_by_d_c0 = edwards97_G1::coeff_d * edwards97_Fq3::non_residue;
    edwards97_twist_mul_by_d_c1 = edwards97_G1::coeff_d;
    edwards97_twist_mul_by_d_c2 = edwards97_G1::coeff_d;
    edwards97_twist_mul_by_q_Y = edwards97_Fq("394610617619779289602567691435");
    edwards97_twist_mul_by_q_Z = edwards97_Fq("394610617619779289602567691435");

    /* choice of group G1 */
    edwards97_G1::G1_zero = edwards97_G1(edwards97_Fq::zero(),
                                     edwards97_Fq::one());
    edwards97_G1::G1_one = edwards97_G1(edwards97_Fq("35853984911660288509491339153"),
                                    edwards97_Fq("383408730566446498428080402051"));
    edwards97_G1::initialized = true;

    // TODO
    edwards97_G1::wnaf_window_table.resize(0);
    edwards97_G1::wnaf_window_table.push_back(9);
    edwards97_G1::wnaf_window_table.push_back(14);
    edwards97_G1::wnaf_window_table.push_back(24);
    edwards97_G1::wnaf_window_table.push_back(117);

    edwards97_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.10]
    edwards97_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.10, 9.69]
    edwards97_G1::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [9.69, 25.21]
    edwards97_G1::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.21, 60.00]
    edwards97_G1::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.00, 149.33]
    edwards97_G1::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [149.33, 369.61]
    edwards97_G1::fixed_base_exp_window_table.push_back(149);
    // window 7 is unbeaten in [369.61, 849.07]
    edwards97_G1::fixed_base_exp_window_table.push_back(370);
    // window 8 is unbeaten in [849.07, 1764.94]
    edwards97_G1::fixed_base_exp_window_table.push_back(849);
    // window 9 is unbeaten in [1764.94, 4429.59]
    edwards97_G1::fixed_base_exp_window_table.push_back(1765);
    // window 10 is unbeaten in [4429.59, 13388.78]
    edwards97_G1::fixed_base_exp_window_table.push_back(4430);
    // window 11 is unbeaten in [13388.78, 15368.00]
    edwards97_G1::fixed_base_exp_window_table.push_back(13389);
    // window 12 is unbeaten in [15368.00, 74912.07]
    edwards97_G1::fixed_base_exp_window_table.push_back(15368);
    // window 13 is unbeaten in [74912.07, 438107.20]
    edwards97_G1::fixed_base_exp_window_table.push_back(74912);
    // window 14 is never the best
    edwards97_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [438107.20, 1045626.18]
    edwards97_G1::fixed_base_exp_window_table.push_back(438107);
    // window 16 is never the best
    edwards97_G1::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [1045626.18, 1577434.48]
    edwards97_G1::fixed_base_exp_window_table.push_back(1045626);
    // window 18 is unbeaten in [1577434.48, 17350594.23]
    edwards97_G1::fixed_base_exp_window_table.push_back(1577434);
    // window 19 is never the best
    edwards97_G1::fixed_base_exp_window_table.push_back(0);
    // window 20 is never the best
    edwards97_G1::fixed_base_exp_window_table.push_back(0);
    // window 21 is unbeaten in [17350594.23, inf]
    edwards97_G1::fixed_base_exp_window_table.push_back(17350594);
    // window 22 is never the best
    edwards97_G1::fixed_base_exp_window_table.push_back(0);

    /* choice of group G2 */

    edwards97_G2::G2_zero = edwards97_G2(edwards97_Fq3::zero(),
                                     edwards97_Fq3::one());
    edwards97_G2::G2_one = edwards97_G2(edwards97_Fq3(edwards97_Fq("123070852602716777894164544077"),
                                                edwards97_Fq("133437722094334994737449128347"),
                                                edwards97_Fq("191829953762958643248750650417")),
                                    edwards97_Fq3(edwards97_Fq("553009370798710752758919914020"),
                                                edwards97_Fq("146270529098674545952140235489"),
                                                edwards97_Fq("251944558233906904530334678700")));
    edwards97_G2::initialized = true;

    // TODO
    edwards97_G2::wnaf_window_table.resize(0);
    edwards97_G2::wnaf_window_table.push_back(6);
    edwards97_G2::wnaf_window_table.push_back(12);
    edwards97_G2::wnaf_window_table.push_back(42);
    edwards97_G2::wnaf_window_table.push_back(97);

    edwards97_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.74]
    edwards97_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.74, 10.67]
    edwards97_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.67, 25.53]
    edwards97_G2::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [25.53, 60.67]
    edwards97_G2::fixed_base_exp_window_table.push_back(26);
    // window 5 is unbeaten in [60.67, 145.77]
    edwards97_G2::fixed_base_exp_window_table.push_back(61);
    // window 6 is unbeaten in [145.77, 356.76]
    edwards97_G2::fixed_base_exp_window_table.push_back(146);
    // window 7 is unbeaten in [356.76, 823.08]
    edwards97_G2::fixed_base_exp_window_table.push_back(357);
    // window 8 is unbeaten in [823.08, 1589.45]
    edwards97_G2::fixed_base_exp_window_table.push_back(823);
    // window 9 is unbeaten in [1589.45, 4135.70]
    edwards97_G2::fixed_base_exp_window_table.push_back(1589);
    // window 10 is unbeaten in [4135.70, 14297.74]
    edwards97_G2::fixed_base_exp_window_table.push_back(4136);
    // window 11 is unbeaten in [14297.74, 16744.85]
    edwards97_G2::fixed_base_exp_window_table.push_back(14298);
    // window 12 is unbeaten in [16744.85, 51768.98]
    edwards97_G2::fixed_base_exp_window_table.push_back(16745);
    // window 13 is unbeaten in [51768.98, 99811.01]
    edwards97_G2::fixed_base_exp_window_table.push_back(51769);
    // window 14 is unbeaten in [99811.01, 193306.72]
    edwards97_G2::fixed_base_exp_window_table.push_back(99811);
    // window 15 is unbeaten in [193306.72, 907184.68]
    edwards97_G2::fixed_base_exp_window_table.push_back(193307);
    // window 16 is never the best
    edwards97_G2::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [907184.68, 1389682.59]
    edwards97_G2::fixed_base_exp_window_table.push_back(907185);
    // window 18 is unbeaten in [1389682.59, 6752695.74]
    edwards97_G2::fixed_base_exp_window_table.push_back(1389683);
    // window 19 is never the best
    edwards97_G2::fixed_base_exp_window_table.push_back(0);
    // window 20 is unbeaten in [6752695.74, 193642894.51]
    edwards97_G2::fixed_base_exp_window_table.push_back(6752696);
    // window 21 is unbeaten in [193642894.51, 226760202.29]
    edwards97_G2::fixed_base_exp_window_table.push_back(193642895);
    // window 22 is unbeaten in [226760202.29, inf]
    edwards97_G2::fixed_base_exp_window_table.push_back(226760202);

    /* pairing parameters */

    edwards97_ate_loop_count = bigint<edwards97_q_limbs>("1356070047940611");
    edwards97_final_exponent_last_chunk_abs_of_w0 = bigint<edwards97_q_limbs>("5424280191762435");
    edwards97_final_exponent_last_chunk_is_w0_neg = true;
    edwards97_final_exponent_last_chunk_w1 = bigint<edwards97_q_limbs>("4");


}
} // libff
