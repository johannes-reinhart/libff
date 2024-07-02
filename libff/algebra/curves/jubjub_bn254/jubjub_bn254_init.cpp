#include <libff/algebra/curves/jubjub_bn254/jubjub_bn254_g1.hpp>
#include <libff/algebra/curves/jubjub_bn254/jubjub_bn254_init.hpp>

namespace libff {

void init_jubjub_bn254_params()
{
    init_jubjub_bn254_fields();

    /* choice of Edwards curve and its twist */

    jubjub_bn254_G1::coeff_a = jubjub_bn254_Fq("168700");
    jubjub_bn254_G1::coeff_d = jubjub_bn254_Fq("168696");

    /* choice of group G1 */
    jubjub_bn254_G1::G1_zero = jubjub_bn254_G1(jubjub_bn254_Fq::zero(),
                                     jubjub_bn254_Fq::one());
    jubjub_bn254_G1::G1_one = jubjub_bn254_G1(jubjub_bn254_Fq("5299619240641551281634865583518297030282874472190772894086521144482721001553"),
                                    jubjub_bn254_Fq("16950150798460657717958625567821834550301663161624707787222815936182638968203"));
    jubjub_bn254_G1::initialized = true;

    // TODO
    jubjub_bn254_G1::wnaf_window_table.resize(0);
    jubjub_bn254_G1::wnaf_window_table.push_back(9);
    jubjub_bn254_G1::wnaf_window_table.push_back(14);
    jubjub_bn254_G1::wnaf_window_table.push_back(24);
    jubjub_bn254_G1::wnaf_window_table.push_back(117);

    jubjub_bn254_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.10]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.10, 9.69]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(4);
    // window 3 is unbeaten in [9.69, 25.21]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.21, 60.00]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [60.00, 149.33]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(60);
    // window 6 is unbeaten in [149.33, 369.61]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(149);
    // window 7 is unbeaten in [369.61, 849.07]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(370);
    // window 8 is unbeaten in [849.07, 1764.94]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(849);
    // window 9 is unbeaten in [1764.94, 4429.59]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(1765);
    // window 10 is unbeaten in [4429.59, 13388.78]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(4430);
    // window 11 is unbeaten in [13388.78, 15368.00]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(13389);
    // window 12 is unbeaten in [15368.00, 74912.07]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(15368);
    // window 13 is unbeaten in [74912.07, 438107.20]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(74912);
    // window 14 is never the best
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [438107.20, 1045626.18]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(438107);
    // window 16 is never the best
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(0);
    // window 17 is unbeaten in [1045626.18, 1577434.48]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(1045626);
    // window 18 is unbeaten in [1577434.48, 17350594.23]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(1577434);
    // window 19 is never the best
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(0);
    // window 20 is never the best
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(0);
    // window 21 is unbeaten in [17350594.23, inf]
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(17350594);
    // window 22 is never the best
    jubjub_bn254_G1::fixed_base_exp_window_table.push_back(0);

}
} // libff
