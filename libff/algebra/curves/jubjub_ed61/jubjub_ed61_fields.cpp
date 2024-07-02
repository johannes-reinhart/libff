#include <libff/algebra/curves/jubjub_ed61/jubjub_ed61_fields.hpp>

namespace libff {

bigint<jubjub_ed61_r_limbs> jubjub_ed61_modulus_r;
bigint<jubjub_ed61_q_limbs> jubjub_ed61_modulus_q;

void init_jubjub_ed61_fields()
{
    using bigint_r = bigint<jubjub_ed61_r_limbs>;
    using bigint_q = bigint<jubjub_ed61_q_limbs>;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    jubjub_ed61_modulus_r = bigint_r("154435754209896193");
    assert(jubjub_ed61_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        jubjub_ed61_Fr::Rsquared = bigint_r("148912062088797258");
        jubjub_ed61_Fr::Rcubed = bigint_r("48817374331921952");
        jubjub_ed61_Fr::inv = 0x2e4dc9fb0cdbe6ff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        jubjub_ed61_Fr::Rsquared = bigint_r("148912062088797258");
        jubjub_ed61_Fr::Rcubed = bigint_r("48817374331921952");
        jubjub_ed61_Fr::inv = 0xcdbe6ff;
    }
    jubjub_ed61_Fr::num_bits = 58;
    jubjub_ed61_Fr::euler = bigint_r("77217877104948096");
    jubjub_ed61_Fr::s = 8;
    jubjub_ed61_Fr::t = bigint_r("603264664882407");
    jubjub_ed61_Fr::t_minus_1_over_2 = bigint_r("301632332441203");
    jubjub_ed61_Fr::multiplicative_generator = jubjub_ed61_Fr("7");
    jubjub_ed61_Fr::root_of_unity = jubjub_ed61_Fr("22108946416905303");
    jubjub_ed61_Fr::nqr = jubjub_ed61_Fr("5");
    jubjub_ed61_Fr::nqr_to_t = jubjub_ed61_Fr("116617487227345212");

    /* parameters for base field Fq */

    jubjub_ed61_modulus_q = bigint_q("1235486033917771777");
    assert(jubjub_ed61_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        jubjub_ed61_Fq::Rsquared = bigint_q("1073255652645550202");
        jubjub_ed61_Fq::Rcubed = bigint_q("357882241080121889");
        jubjub_ed61_Fq::inv = 0x5084f000809fffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        jubjub_ed61_Fq::Rsquared = bigint_q("1073255652645550202");
        jubjub_ed61_Fq::Rcubed = bigint_q("357882241080121889");
        jubjub_ed61_Fq::inv = 0x809fffff;
    }
    jubjub_ed61_Fq::num_bits = 61;
    jubjub_ed61_Fq::euler = bigint_q("617743016958885888");
    jubjub_ed61_Fq::s = 21;
    jubjub_ed61_Fq::t = bigint_q("589125649413");
    jubjub_ed61_Fq::t_minus_1_over_2 = bigint_q("294562824706");
    jubjub_ed61_Fq::multiplicative_generator = jubjub_ed61_Fq("13");
    jubjub_ed61_Fq::root_of_unity = jubjub_ed61_Fq("909855932006627559");
    jubjub_ed61_Fq::nqr = jubjub_ed61_Fq("5");
    jubjub_ed61_Fq::nqr_to_t = jubjub_ed61_Fq("468172926162659554");

    /* parameters for twist field Fq3 */

    jubjub_ed61_Fq3::euler = bigint<3*jubjub_ed61_q_limbs>("942938841794923322102999583456121473051897856672137216");
    jubjub_ed61_Fq3::s = 21;
    jubjub_ed61_Fq3::t = bigint<3*jubjub_ed61_q_limbs>("899256555361674615958213408905145142604730469391");
    jubjub_ed61_Fq3::t_minus_1_over_2 = bigint<3*jubjub_ed61_q_limbs>("449628277680837307979106704452572571302365234695");
    jubjub_ed61_Fq3::non_residue = jubjub_ed61_Fq("13");
    jubjub_ed61_Fq3::nqr = jubjub_ed61_Fq3(jubjub_ed61_Fq("5"),jubjub_ed61_Fq("0"),jubjub_ed61_Fq("0"));
    jubjub_ed61_Fq3::nqr_to_t = jubjub_ed61_Fq3(jubjub_ed61_Fq("743359549615257709"),jubjub_ed61_Fq("0"),jubjub_ed61_Fq("0"));
    jubjub_ed61_Fq3::Frobenius_coeffs_c1[0] = jubjub_ed61_Fq("1");
    jubjub_ed61_Fq3::Frobenius_coeffs_c1[1] = jubjub_ed61_Fq("4007657475");
    jubjub_ed61_Fq3::Frobenius_coeffs_c1[2] = jubjub_ed61_Fq("1235486029910114301");
    jubjub_ed61_Fq3::Frobenius_coeffs_c2[0] = jubjub_ed61_Fq("1");
    jubjub_ed61_Fq3::Frobenius_coeffs_c2[1] = jubjub_ed61_Fq("1235486029910114301");
    jubjub_ed61_Fq3::Frobenius_coeffs_c2[2] = jubjub_ed61_Fq("4007657475");

    /* parameters for Fq6 */

    jubjub_ed61_Fq6::euler = bigint<6*jubjub_ed61_q_limbs>("1778267318731102867374494615436930276822520310043298290702743492661414821401862074046723012800168171608735744");
    jubjub_ed61_Fq6::s = 22;
    jubjub_ed61_Fq6::t = bigint<6*jubjub_ed61_q_limbs>("847943934789229806601760204046692980204830317517899651862498995142657671643191372893678194427570424847");
    jubjub_ed61_Fq6::t_minus_1_over_2 = bigint<6*jubjub_ed61_q_limbs>("423971967394614903300880102023346490102415158758949825931249497571328835821595686446839097213785212423");
    jubjub_ed61_Fq6::non_residue = jubjub_ed61_Fq("13");
    jubjub_ed61_Fq6::nqr = jubjub_ed61_Fq6(jubjub_ed61_Fq3(jubjub_ed61_Fq("3"),jubjub_ed61_Fq("0"),jubjub_ed61_Fq("0")),jubjub_ed61_Fq3::one());
    jubjub_ed61_Fq6::nqr_to_t = jubjub_ed61_Fq6(jubjub_ed61_Fq3::zero(),jubjub_ed61_Fq3(jubjub_ed61_Fq("0"),jubjub_ed61_Fq("1067377170496400159"),jubjub_ed61_Fq("0")));
    jubjub_ed61_Fq6::Frobenius_coeffs_c1[0] = jubjub_ed61_Fq("1");
    jubjub_ed61_Fq6::Frobenius_coeffs_c1[1] = jubjub_ed61_Fq("4007657476");
    jubjub_ed61_Fq6::Frobenius_coeffs_c1[2] = jubjub_ed61_Fq("4007657475");
    jubjub_ed61_Fq6::Frobenius_coeffs_c1[3] = jubjub_ed61_Fq("1235486033917771776");
    jubjub_ed61_Fq6::Frobenius_coeffs_c1[4] = jubjub_ed61_Fq("1235486029910114301");
    jubjub_ed61_Fq6::Frobenius_coeffs_c1[5] = jubjub_ed61_Fq("1235486029910114302");
    jubjub_ed61_Fq6::my_Fp2::non_residue = jubjub_ed61_Fq3::non_residue;
}

} // namespace libff
