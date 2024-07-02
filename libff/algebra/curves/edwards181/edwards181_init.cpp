#include <libff/algebra/curves/edwards181/edwards181_g1.hpp>
#include <libff/algebra/curves/edwards181/edwards181_g2.hpp>
#include <libff/algebra/curves/edwards181/edwards181_init.hpp>

namespace libff {

edwards181_Fq3 edwards181_twist;
edwards181_Fq3 edwards181_twist_coeff_a;
edwards181_Fq3 edwards181_twist_coeff_d;
edwards181_Fq edwards181_twist_mul_by_a_c0;
edwards181_Fq edwards181_twist_mul_by_a_c1;
edwards181_Fq edwards181_twist_mul_by_a_c2;
edwards181_Fq edwards181_twist_mul_by_d_c0;
edwards181_Fq edwards181_twist_mul_by_d_c1;
edwards181_Fq edwards181_twist_mul_by_d_c2;
edwards181_Fq edwards181_twist_mul_by_q_Y;
edwards181_Fq edwards181_twist_mul_by_q_Z;

bigint<edwards181_q_limbs> edwards181_ate_loop_count;
bigint<edwards181_q_limbs> edwards181_final_exponent_last_chunk_abs_of_w0;
bool edwards181_final_exponent_last_chunk_is_w0_neg;
bigint<edwards181_q_limbs> edwards181_final_exponent_last_chunk_w1;

void init_edwards181_params()
{
    init_edwards181_fields();

    /* choice of Edwards curve and its twist */

    edwards181_G1::coeff_a = edwards181_Fq("1");
    edwards181_G1::coeff_d = edwards181_Fq("600581931845324488256649384912508268813600056237543024");
    edwards181_twist = edwards181_Fq3(edwards181_Fq::zero(), edwards181_Fq::one(), edwards181_Fq::zero());
    edwards181_twist_coeff_a = edwards181_G1::coeff_a * edwards181_twist;
    edwards181_twist_coeff_d = edwards181_G1::coeff_d * edwards181_twist;
    edwards181_twist_mul_by_a_c0 = edwards181_G1::coeff_a * edwards181_Fq3::non_residue;
    edwards181_twist_mul_by_a_c1 = edwards181_G1::coeff_a;
    edwards181_twist_mul_by_a_c2 = edwards181_G1::coeff_a;
    edwards181_twist_mul_by_d_c0 = edwards181_G1::coeff_d * edwards181_Fq3::non_residue;
    edwards181_twist_mul_by_d_c1 = edwards181_G1::coeff_d;
    edwards181_twist_mul_by_d_c2 = edwards181_G1::coeff_d;
    edwards181_twist_mul_by_q_Y = edwards181_Fq("1073752683758513276629212192812154536507607213288832062");
    edwards181_twist_mul_by_q_Z = edwards181_Fq("1073752683758513276629212192812154536507607213288832062");

    /* choice of group G1 */
    edwards181_G1::G1_zero = edwards181_G1(edwards181_Fq::zero(),
                                     edwards181_Fq::one());
    edwards181_G1::G1_one = edwards181_G1(edwards181_Fq("3713709671941291996998665608188072510389821008693530490"),
                                    edwards181_Fq("4869953702976555123067178261685365085639705297852816679"));
    edwards181_G1::initialized = true;

    // TODO
    edwards181_G1::wnaf_window_table.resize(0);
    edwards181_G1::wnaf_window_table.push_back(9);
    edwards181_G1::wnaf_window_table.push_back(14);
    edwards181_G1::wnaf_window_table.push_back(24);
    edwards181_G1::wnaf_window_table.push_back(117);

    edwards181_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.10]
    edwards181_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.10, 9.69]
    edwards181_G1::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [9.69, 25.21]
    edwards181_G1::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.21, 60.00]
    edwards181_G1::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.00, 149.33]
    edwards181_G1::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [149.33, 369.61]
    edwards181_G1::fixed_base_exp_window_table.push_back(149);
    // window 7 is unbeaten in [369.61, 849.07]
    edwards181_G1::fixed_base_exp_window_table.push_back(370);
    // window 8 is unbeaten in [849.07, 1764.94]
    edwards181_G1::fixed_base_exp_window_table.push_back(849);
    // window 9 is unbeaten in [1764.94, 4429.59]
    edwards181_G1::fixed_base_exp_window_table.push_back(1765);
    // window 10 is unbeaten in [4429.59, 13388.78]
    edwards181_G1::fixed_base_exp_window_table.push_back(4430);
    // window 11 is unbeaten in [13388.78, 15368.00]
    edwards181_G1::fixed_base_exp_window_table.push_back(13389);
    // window 12 is unbeaten in [15368.00, 74912.07]
    edwards181_G1::fixed_base_exp_window_table.push_back(15368);
    // window 13 is unbeaten in [74912.07, 438107.20]
    edwards181_G1::fixed_base_exp_window_table.push_back(74912);
    // window 14 is never the best
    edwards181_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [438107.20, 1045626.18]
    edwards181_G1::fixed_base_exp_window_table.push_back(438107);
    // window 16 is never the best
    edwards181_G1::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [1045626.18, 1577434.48]
    edwards181_G1::fixed_base_exp_window_table.push_back(1045626);
    // window 18 is unbeaten in [1577434.48, 17350594.23]
    edwards181_G1::fixed_base_exp_window_table.push_back(1577434);
    // window 19 is never the best
    edwards181_G1::fixed_base_exp_window_table.push_back(0);
    // window 20 is never the best
    edwards181_G1::fixed_base_exp_window_table.push_back(0);
    // window 21 is unbeaten in [17350594.23, inf]
    edwards181_G1::fixed_base_exp_window_table.push_back(17350594);
    // window 22 is never the best
    edwards181_G1::fixed_base_exp_window_table.push_back(0);

    /* choice of group G2 */

    edwards181_G2::G2_zero = edwards181_G2(edwards181_Fq3::zero(),
                                     edwards181_Fq3::one());
    edwards181_G2::G2_one = edwards181_G2(edwards181_Fq3(edwards181_Fq("2985169216050045936447045373806587563703160287265791273"),
                                                edwards181_Fq("4056487470863548806908274038624319453234070766729386632"),
                                                edwards181_Fq("5420687431981356209724320472992390796720155111208612181")),
                                    edwards181_Fq3(edwards181_Fq("5511769377402467244479502559508855524462271277798865946"),
                                                edwards181_Fq("4596331147700187028340200824898965930392167581458537157"),
                                                edwards181_Fq("2612300818842229264542909713627888465410600578484324309")));
    edwards181_G2::initialized = true;

    // TODO
    edwards181_G2::wnaf_window_table.resize(0);
    edwards181_G2::wnaf_window_table.push_back(6);
    edwards181_G2::wnaf_window_table.push_back(12);
    edwards181_G2::wnaf_window_table.push_back(42);
    edwards181_G2::wnaf_window_table.push_back(97);

    edwards181_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.74]
    edwards181_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.74, 10.67]
    edwards181_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.67, 25.53]
    edwards181_G2::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [25.53, 60.67]
    edwards181_G2::fixed_base_exp_window_table.push_back(26);
    // window 5 is unbeaten in [60.67, 145.77]
    edwards181_G2::fixed_base_exp_window_table.push_back(61);
    // window 6 is unbeaten in [145.77, 356.76]
    edwards181_G2::fixed_base_exp_window_table.push_back(146);
    // window 7 is unbeaten in [356.76, 823.08]
    edwards181_G2::fixed_base_exp_window_table.push_back(357);
    // window 8 is unbeaten in [823.08, 1589.45]
    edwards181_G2::fixed_base_exp_window_table.push_back(823);
    // window 9 is unbeaten in [1589.45, 4135.70]
    edwards181_G2::fixed_base_exp_window_table.push_back(1589);
    // window 10 is unbeaten in [4135.70, 14297.74]
    edwards181_G2::fixed_base_exp_window_table.push_back(4136);
    // window 11 is unbeaten in [14297.74, 16744.85]
    edwards181_G2::fixed_base_exp_window_table.push_back(14298);
    // window 12 is unbeaten in [16744.85, 51768.98]
    edwards181_G2::fixed_base_exp_window_table.push_back(16745);
    // window 13 is unbeaten in [51768.98, 99811.01]
    edwards181_G2::fixed_base_exp_window_table.push_back(51769);
    // window 14 is unbeaten in [99811.01, 193306.72]
    edwards181_G2::fixed_base_exp_window_table.push_back(99811);
    // window 15 is unbeaten in [193306.72, 907184.68]
    edwards181_G2::fixed_base_exp_window_table.push_back(193307);
    // window 16 is never the best
    edwards181_G2::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [907184.68, 1389682.59]
    edwards181_G2::fixed_base_exp_window_table.push_back(907185);
    // window 18 is unbeaten in [1389682.59, 6752695.74]
    edwards181_G2::fixed_base_exp_window_table.push_back(1389683);
    // window 19 is never the best
    edwards181_G2::fixed_base_exp_window_table.push_back(0);
    // window 20 is unbeaten in [6752695.74, 193642894.51]
    edwards181_G2::fixed_base_exp_window_table.push_back(6752696);
    // window 21 is unbeaten in [193642894.51, 226760202.29]
    edwards181_G2::fixed_base_exp_window_table.push_back(193642895);
    // window 22 is unbeaten in [226760202.29, inf]
    edwards181_G2::fixed_base_exp_window_table.push_back(226760202);

    /* pairing parameters */

    edwards181_ate_loop_count = bigint<edwards181_q_limbs>("4492509698523932320491110403");
    edwards181_final_exponent_last_chunk_abs_of_w0 = bigint<edwards181_q_limbs>("17970038794095729281964441603");
    edwards181_final_exponent_last_chunk_is_w0_neg = true;
    edwards181_final_exponent_last_chunk_w1 = bigint<edwards181_q_limbs>("4");


}
} // libff
