#include <libff/algebra/curves/jubjub_bn124/jubjub_bn124_fields.hpp>

namespace libff {

bigint<jubjub_bn124_r_limbs> jubjub_bn124_modulus_r;
bigint<jubjub_bn124_q_limbs> jubjub_bn124_modulus_q;

void init_jubjub_bn124_fields()
{
    using bigint_r = bigint<jubjub_bn124_r_limbs>;
    using bigint_q = bigint<jubjub_bn124_q_limbs>;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    jubjub_bn124_modulus_r = bigint_r("2125016665599104007263567884627882177");
    assert(jubjub_bn124_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        jubjub_bn124_Fr::Rsquared = bigint_r("2118178146368725646386118482215374336");
        jubjub_bn124_Fr::Rcubed = bigint_r("1185055689807017249812904351266995641");
        jubjub_bn124_Fr::inv = 0x61d1fac12f51a0bf;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        jubjub_bn124_Fr::Rsquared = bigint_r("2118178146368725646386118482215374336");
        jubjub_bn124_Fr::Rcubed = bigint_r("1185055689807017249812904351266995641");
        jubjub_bn124_Fr::inv = 0x2f51a0bf;
    }
    jubjub_bn124_Fr::num_bits = 121;
    jubjub_bn124_Fr::euler = bigint_r("1062508332799552003631783942313941088");
    jubjub_bn124_Fr::s = 6;
    jubjub_bn124_Fr::t = bigint_r("33203385399986000113493248197310659");
    jubjub_bn124_Fr::t_minus_1_over_2 = bigint_r("16601692699993000056746624098655329");
    jubjub_bn124_Fr::multiplicative_generator = jubjub_bn124_Fr("5");
    jubjub_bn124_Fr::root_of_unity = jubjub_bn124_Fr("1054393897166881011445462479212472982");
    jubjub_bn124_Fr::nqr = jubjub_bn124_Fr("5");
    jubjub_bn124_Fr::nqr_to_t = jubjub_bn124_Fr("1054393897166881011445462479212472982");

    /* parameters for base field Fq */

    jubjub_bn124_modulus_q = bigint_q("17000133324792832058895897937997463553");
    assert(jubjub_bn124_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        jubjub_bn124_Fq::Rsquared = bigint_q("11418104352409030789881016589005285164");
        jubjub_bn124_Fq::Rcubed = bigint_q("14736843594676410769816466500475533132");
        jubjub_bn124_Fq::inv = 0x6f1a7d37f9ffffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        jubjub_bn124_Fq::Rsquared = bigint_q("11418104352409030789881016589005285164");
        jubjub_bn124_Fq::Rcubed = bigint_q("14736843594676410769816466500475533132");
        jubjub_bn124_Fq::inv = 0xf9ffffff;
    }
    jubjub_bn124_Fq::num_bits = 124;
    jubjub_bn124_Fq::euler = bigint_q("8500066662396416029447948968998731776");
    jubjub_bn124_Fq::s = 25;
    jubjub_bn124_Fq::t = bigint_q("506643453979278566208359537661");
    jubjub_bn124_Fq::t_minus_1_over_2 = bigint_q("253321726989639283104179768830");
    jubjub_bn124_Fq::multiplicative_generator = jubjub_bn124_Fq("7");
    jubjub_bn124_Fq::root_of_unity = jubjub_bn124_Fq("7041117609370421414184985844533865071");
    jubjub_bn124_Fq::nqr = jubjub_bn124_Fq("5");
    jubjub_bn124_Fq::nqr_to_t = jubjub_bn124_Fq("5402570782257262925316013050580052157");

    /* parameters for twist field Fq3 */

    jubjub_bn124_Fq3::euler = bigint<3*jubjub_bn124_q_limbs>("2456557796750969142273466118625016884658446119321628675722106284372727175931587668035076640293019332490913906688");
    jubjub_bn124_Fq3::s = 25;
    jubjub_bn124_Fq3::t = bigint<3*jubjub_bn124_q_limbs>("146422254845557757751552231229842715541031725366212646706229822896285484786724309208099641817392071037943");
    jubjub_bn124_Fq3::t_minus_1_over_2 = bigint<3*jubjub_bn124_q_limbs>("73211127422778878875776115614921357770515862683106323353114911448142742393362154604049820908696035518971");
    jubjub_bn124_Fq3::non_residue = jubjub_bn124_Fq("7");
    jubjub_bn124_Fq3::nqr = jubjub_bn124_Fq3(jubjub_bn124_Fq("5"),jubjub_bn124_Fq("0"),jubjub_bn124_Fq("0"));
    jubjub_bn124_Fq3::nqr_to_t = jubjub_bn124_Fq3(jubjub_bn124_Fq("13346094514203034772038987226297172149"),jubjub_bn124_Fq("0"),jubjub_bn124_Fq("0"));
    jubjub_bn124_Fq3::Frobenius_coeffs_c1[0] = jubjub_bn124_Fq("1");
    jubjub_bn124_Fq3::Frobenius_coeffs_c1[1] = jubjub_bn124_Fq("17000133304285230517543828587634978595");
    jubjub_bn124_Fq3::Frobenius_coeffs_c1[2] = jubjub_bn124_Fq("20507601541352069350362484957");
    jubjub_bn124_Fq3::Frobenius_coeffs_c2[0] = jubjub_bn124_Fq("1");
    jubjub_bn124_Fq3::Frobenius_coeffs_c2[1] = jubjub_bn124_Fq("20507601541352069350362484957");
    jubjub_bn124_Fq3::Frobenius_coeffs_c2[2] = jubjub_bn124_Fq("17000133304285230517543828587634978595");

    /* parameters for Fq6 */

    jubjub_bn124_Fq6::euler = bigint<6*jubjub_bn124_q_limbs>("12069352417555951637157174735894474908586518481163501213479245201132670232033110524716009133688055005142313160411064967370888635034976948596897759989138505065578172849773744493873044759854289602465935071491335888330570072064");
    jubjub_bn124_Fq6::s = 26;
    jubjub_bn124_Fq6::t = bigint<6*jubjub_bn124_q_limbs>("359694731758712280903970442291929570096329405342444813653208172355075783492121414086699757983924597655007635367246418218937177510111837047246031760845735820102041150622777476724178932900854635312138052925209280500727");
    jubjub_bn124_Fq6::t_minus_1_over_2 = bigint<6*jubjub_bn124_q_limbs>("179847365879356140451985221145964785048164702671222406826604086177537891746060707043349878991962298827503817683623209109468588755055918523623015880422867910051020575311388738362089466450427317656069026462604640250363");
    jubjub_bn124_Fq6::non_residue = jubjub_bn124_Fq("7");
    jubjub_bn124_Fq6::nqr = jubjub_bn124_Fq6(jubjub_bn124_Fq3(jubjub_bn124_Fq("2"),jubjub_bn124_Fq("0"),jubjub_bn124_Fq("0")),jubjub_bn124_Fq3::one());
    jubjub_bn124_Fq6::nqr_to_t = jubjub_bn124_Fq6(jubjub_bn124_Fq3::zero(),jubjub_bn124_Fq3(jubjub_bn124_Fq("0"),jubjub_bn124_Fq("13038214925745496114411627393780209517"),jubjub_bn124_Fq("0")));
    jubjub_bn124_Fq6::Frobenius_coeffs_c1[0] = jubjub_bn124_Fq("1");
    jubjub_bn124_Fq6::Frobenius_coeffs_c1[1] = jubjub_bn124_Fq("17000133304285230517543828587634978596");
    jubjub_bn124_Fq6::Frobenius_coeffs_c1[2] = jubjub_bn124_Fq("17000133304285230517543828587634978595");
    jubjub_bn124_Fq6::Frobenius_coeffs_c1[3] = jubjub_bn124_Fq("17000133324792832058895897937997463552");
    jubjub_bn124_Fq6::Frobenius_coeffs_c1[4] = jubjub_bn124_Fq("20507601541352069350362484957");
    jubjub_bn124_Fq6::Frobenius_coeffs_c1[5] = jubjub_bn124_Fq("20507601541352069350362484958");
    jubjub_bn124_Fq6::my_Fp2::non_residue = jubjub_bn124_Fq3::non_residue;
}

} // namespace libff
