#include <libff/algebra/curves/jubjub_ed58/jubjub_ed58_g1.hpp>
#include <libff/algebra/curves/jubjub_ed58/jubjub_ed58_init.hpp>

namespace libff {

void init_jubjub_ed58_params()
{
    init_jubjub_ed58_fields();

    /* choice of Edwards curve and its twist */

    jubjub_ed58_G1::coeff_a = jubjub_ed58_Fq("78956");
    jubjub_ed58_G1::coeff_d = jubjub_ed58_Fq("78952");

    /* choice of group G1 */
    jubjub_ed58_G1::G1_zero = jubjub_ed58_G1(jubjub_ed58_Fq::zero(),
                                     jubjub_ed58_Fq::one());
    jubjub_ed58_G1::G1_one = jubjub_ed58_G1(jubjub_ed58_Fq("136172154314006764"),
                                    jubjub_ed58_Fq("153203390242712599"));
    jubjub_ed58_G1::initialized = true;

    // TODO
    jubjub_ed58_G1::wnaf_window_table.resize(0);
    jubjub_ed58_G1::wnaf_window_table.push_back(9);
    jubjub_ed58_G1::wnaf_window_table.push_back(14);
    jubjub_ed58_G1::wnaf_window_table.push_back(24);
    jubjub_ed58_G1::wnaf_window_table.push_back(117);

    jubjub_ed58_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.10]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.10, 9.69]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [9.69, 25.21]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.21, 60.00]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.00, 149.33]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [149.33, 369.61]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(149);
    // window 7 is unbeaten in [369.61, 849.07]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(370);
    // window 8 is unbeaten in [849.07, 1764.94]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(849);
    // window 9 is unbeaten in [1764.94, 4429.59]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(1765);
    // window 10 is unbeaten in [4429.59, 13388.78]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(4430);
    // window 11 is unbeaten in [13388.78, 15368.00]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(13389);
    // window 12 is unbeaten in [15368.00, 74912.07]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(15368);
    // window 13 is unbeaten in [74912.07, 438107.20]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(74912);
    // window 14 is never the best
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [438107.20, 1045626.18]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(438107);
    // window 16 is never the best
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [1045626.18, 1577434.48]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(1045626);
    // window 18 is unbeaten in [1577434.48, 17350594.23]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(1577434);
    // window 19 is never the best
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(0);
    // window 20 is never the best
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(0);
    // window 21 is unbeaten in [17350594.23, inf]
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(17350594);
    // window 22 is never the best
    jubjub_ed58_G1::fixed_base_exp_window_table.push_back(0);

}
} // libff
