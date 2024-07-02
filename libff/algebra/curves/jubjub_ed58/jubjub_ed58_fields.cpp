#include <libff/algebra/curves/jubjub_ed58/jubjub_ed58_fields.hpp>

namespace libff {

bigint<jubjub_ed58_r_limbs> jubjub_ed58_modulus_r;
bigint<jubjub_ed58_q_limbs> jubjub_ed58_modulus_q;

void init_jubjub_ed58_fields()
{
    using bigint_r = bigint<jubjub_ed58_r_limbs>;
    using bigint_q = bigint<jubjub_ed58_q_limbs>;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    jubjub_ed58_modulus_r = bigint_r("26375806633482721");
    assert(jubjub_ed58_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        jubjub_ed58_Fr::Rsquared = bigint_r("14579298408581276");
        jubjub_ed58_Fr::Rcubed = bigint_r("2435545966331503");
        jubjub_ed58_Fr::inv = 0x2241eeec0b7795df;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        jubjub_ed58_Fr::Rsquared = bigint_r("14579298408581276");
        jubjub_ed58_Fr::Rcubed = bigint_r("2435545966331503");
        jubjub_ed58_Fr::inv = 0xb7795df;
    }
    jubjub_ed58_Fr::num_bits = 55;
    jubjub_ed58_Fr::euler = bigint_r("13187903316741360");
    jubjub_ed58_Fr::s = 5;
    jubjub_ed58_Fr::t = bigint_r("824243957296335");
    jubjub_ed58_Fr::t_minus_1_over_2 = bigint_r("412121978648167");
    jubjub_ed58_Fr::multiplicative_generator = jubjub_ed58_Fr("22");
    jubjub_ed58_Fr::root_of_unity = jubjub_ed58_Fr("22336128816477168");
    jubjub_ed58_Fr::nqr = jubjub_ed58_Fr("11");
    jubjub_ed58_Fr::nqr_to_t = jubjub_ed58_Fr("11568150827995564");

    /* parameters for base field Fq */

    jubjub_ed58_modulus_q = bigint_q("211006452744585217");
    assert(jubjub_ed58_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        jubjub_ed58_Fq::Rsquared = bigint_q("69518267247397641");
        jubjub_ed58_Fq::Rcubed = bigint_q("193528881679320662");
        jubjub_ed58_Fq::inv = 0x8de40f003527ffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        jubjub_ed58_Fq::Rsquared = bigint_q("69518267247397641");
        jubjub_ed58_Fq::Rcubed = bigint_q("193528881679320662");
        jubjub_ed58_Fq::inv = 0x3527ffff;
    }
    jubjub_ed58_Fq::num_bits = 58;
    jubjub_ed58_Fq::euler = bigint_q("105503226372292608");
    jubjub_ed58_Fq::s = 19;
    jubjub_ed58_Fq::t = bigint_q("402462869157");
    jubjub_ed58_Fq::t_minus_1_over_2 = bigint_q("201231434578");
    jubjub_ed58_Fq::multiplicative_generator = jubjub_ed58_Fq("7");
    jubjub_ed58_Fq::root_of_unity = jubjub_ed58_Fq("101722017355283520");
    jubjub_ed58_Fq::nqr = jubjub_ed58_Fq("5");
    jubjub_ed58_Fq::nqr_to_t = jubjub_ed58_Fq("85346003589280634");

    /* parameters for twist field Fq3 */

    jubjub_ed58_Fq3::euler = bigint<3*jubjub_ed58_q_limbs>("4697396437141051372361274092695161632321381387206656");
    jubjub_ed58_Fq3::s = 19;
    jubjub_ed58_Fq3::t = bigint<3*jubjub_ed58_q_limbs>("17919145344318585862584205981045385865483785199");
    jubjub_ed58_Fq3::t_minus_1_over_2 = bigint<3*jubjub_ed58_q_limbs>("8959572672159292931292102990522692932741892599");
    jubjub_ed58_Fq3::non_residue = jubjub_ed58_Fq("7");
    jubjub_ed58_Fq3::nqr = jubjub_ed58_Fq3(jubjub_ed58_Fq("5"),jubjub_ed58_Fq("0"),jubjub_ed58_Fq("0"));
    jubjub_ed58_Fq3::nqr_to_t = jubjub_ed58_Fq3(jubjub_ed58_Fq("105725275087470683"),jubjub_ed58_Fq("0"),jubjub_ed58_Fq("0"));
    jubjub_ed58_Fq3::Frobenius_coeffs_c1[0] = jubjub_ed58_Fq("1");
    jubjub_ed58_Fq3::Frobenius_coeffs_c1[1] = jubjub_ed58_Fq("1656225795");
    jubjub_ed58_Fq3::Frobenius_coeffs_c1[2] = jubjub_ed58_Fq("211006451088359421");
    jubjub_ed58_Fq3::Frobenius_coeffs_c2[0] = jubjub_ed58_Fq("1");
    jubjub_ed58_Fq3::Frobenius_coeffs_c2[1] = jubjub_ed58_Fq("211006451088359421");
    jubjub_ed58_Fq3::Frobenius_coeffs_c2[2] = jubjub_ed58_Fq("1656225795");

    /* parameters for Fq6 */

    jubjub_ed58_Fq6::euler = bigint<6*jubjub_ed58_q_limbs>("44131066575330886793895027778516637509254112712838848540541689628929467860280177007299381748495675817984");
    jubjub_ed58_Fq6::s = 20;
    jubjub_ed58_Fq6::t = bigint<6*jubjub_ed58_q_limbs>("84173329497014783466138892704995417612560487199475953179438952691897330971298555388067973610869743");
    jubjub_ed58_Fq6::t_minus_1_over_2 = bigint<6*jubjub_ed58_q_limbs>("42086664748507391733069446352497708806280243599737976589719476345948665485649277694033986805434871");
    jubjub_ed58_Fq6::non_residue = jubjub_ed58_Fq("7");
    jubjub_ed58_Fq6::nqr = jubjub_ed58_Fq6(jubjub_ed58_Fq3(jubjub_ed58_Fq("4"),jubjub_ed58_Fq("0"),jubjub_ed58_Fq("0")),jubjub_ed58_Fq3::one());
    jubjub_ed58_Fq6::nqr_to_t = jubjub_ed58_Fq6(jubjub_ed58_Fq3::zero(),jubjub_ed58_Fq3(jubjub_ed58_Fq("0"),jubjub_ed58_Fq("89591299763421710"),jubjub_ed58_Fq("0")));
    jubjub_ed58_Fq6::Frobenius_coeffs_c1[0] = jubjub_ed58_Fq("1");
    jubjub_ed58_Fq6::Frobenius_coeffs_c1[1] = jubjub_ed58_Fq("1656225796");
    jubjub_ed58_Fq6::Frobenius_coeffs_c1[2] = jubjub_ed58_Fq("1656225795");
    jubjub_ed58_Fq6::Frobenius_coeffs_c1[3] = jubjub_ed58_Fq("211006452744585216");
    jubjub_ed58_Fq6::Frobenius_coeffs_c1[4] = jubjub_ed58_Fq("211006451088359421");
    jubjub_ed58_Fq6::Frobenius_coeffs_c1[5] = jubjub_ed58_Fq("211006451088359422");
    jubjub_ed58_Fq6::my_Fp2::non_residue = jubjub_ed58_Fq3::non_residue;
}

} // namespace libff
