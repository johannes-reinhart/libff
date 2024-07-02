#include <libff/algebra/curves/edwards97/edwards97_fields.hpp>

namespace libff {

bigint<edwards97_r_limbs> edwards97_modulus_r;
bigint<edwards97_q_limbs> edwards97_modulus_q;

void init_edwards97_fields()
{
    using bigint_r = bigint<edwards97_r_limbs>;
    using bigint_q = bigint<edwards97_q_limbs>;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    edwards97_modulus_r = bigint_r("141455844224742490147094691841");
    assert(edwards97_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        edwards97_Fr::Rsquared = bigint_r("75217608311129884939407662890");
        edwards97_Fr::Rcubed = bigint_r("128645697791147296497719755611");
        edwards97_Fr::inv = 0x69fd1db9f8ce7fff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        edwards97_Fr::Rsquared = bigint_r("75217608311129884939407662890");
        edwards97_Fr::Rcubed = bigint_r("128645697791147296497719755611");
        edwards97_Fr::inv = 0xf8ce7fff;
    }
    edwards97_Fr::num_bits = 97;
    edwards97_Fr::euler = bigint_r("70727922112371245073547345920");
    edwards97_Fr::s = 15;
    edwards97_Fr::t = bigint_r("4316889777366409001071005");
    edwards97_Fr::t_minus_1_over_2 = bigint_r("2158444888683204500535502");
    edwards97_Fr::multiplicative_generator = edwards97_Fr("13");
    edwards97_Fr::root_of_unity = edwards97_Fr("111466351910010104363354434347");
    edwards97_Fr::nqr = edwards97_Fr("13");
    edwards97_Fr::nqr_to_t = edwards97_Fr("111466351910010104363354434347");

    /* parameters for base field Fq */

    edwards97_modulus_q = bigint_q("565823376898968604518330826753");
    assert(edwards97_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        edwards97_Fq::Rsquared = bigint_q("526869161015607510292097823780");
        edwards97_Fq::Rcubed = bigint_q("496700202014202979788481684526");
        edwards97_Fq::inv = 0xaf6c2afbf9ba7fff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        edwards97_Fq::Rsquared = bigint_q("526869161015607510292097823780");
        edwards97_Fq::Rcubed = bigint_q("496700202014202979788481684526");
        edwards97_Fq::inv = 0xf9ba7fff;
    }
    edwards97_Fq::num_bits = 99;
    edwards97_Fq::euler = bigint_q("282911688449484302259165413376");
    edwards97_Fq::s = 15;
    edwards97_Fq::t = bigint_q("17267559109465594620310389");
    edwards97_Fq::t_minus_1_over_2 = bigint_q("8633779554732797310155194");
    edwards97_Fq::multiplicative_generator = edwards97_Fq("5");
    edwards97_Fq::root_of_unity = edwards97_Fq("404024150737442750347811534883");
    edwards97_Fq::nqr = edwards97_Fq("5");
    edwards97_Fq::nqr_to_t = edwards97_Fq("404024150737442750347811534883");

    /* parameters for twist field Fq3 */

    edwards97_Fq3::euler = bigint<3*edwards97_q_limbs>("90575901077180162453644222721688651680972959905418818484112410359588421548292203934629888");
    edwards97_Fq3::s = 15;
    edwards97_Fq3::t = bigint<3*edwards97_q_limbs>("5528314274730234524758558515728067119200009759852222807868189108861597994890881587807");
    edwards97_Fq3::t_minus_1_over_2 = bigint<3*edwards97_q_limbs>("2764157137365117262379279257864033559600004879926111403934094554430798997445440793903");
    edwards97_Fq3::non_residue = edwards97_Fq("5");
    edwards97_Fq3::nqr = edwards97_Fq3(edwards97_Fq("5"),edwards97_Fq("0"),edwards97_Fq("0"));
    edwards97_Fq3::nqr_to_t = edwards97_Fq3(edwards97_Fq("179411943930830432587917996150"),edwards97_Fq("0"),edwards97_Fq("0"));
    edwards97_Fq3::Frobenius_coeffs_c1[0] = edwards97_Fq("1");
    edwards97_Fq3::Frobenius_coeffs_c1[1] = edwards97_Fq("394610617619779289602567691434");
    edwards97_Fq3::Frobenius_coeffs_c1[2] = edwards97_Fq("171212759279189314915763135318");
    edwards97_Fq3::Frobenius_coeffs_c2[0] = edwards97_Fq("1");
    edwards97_Fq3::Frobenius_coeffs_c2[1] = edwards97_Fq("171212759279189314915763135318");
    edwards97_Fq3::Frobenius_coeffs_c2[2] = edwards97_Fq("394610617619779289602567691434");

    /* parameters for Fq6 */

    edwards97_Fq6::euler = bigint<6*edwards97_q_limbs>("16407987711886253026174339778588314281599367419007115232539561695740466230132674128732333363932787412972054263604013850669983501136804399772397158168017091471775968303518955044864");
    edwards97_Fq6::s = 16;
    edwards97_Fq6::t = bigint<6*edwards97_q_limbs>("500732046871528717839793084063364083300761945160129249039903616203017157901998111838755290647362897124391304431274836751403305088403454582897862492920443465325194345200163423");
    edwards97_Fq6::t_minus_1_over_2 = bigint<6*edwards97_q_limbs>("250366023435764358919896542031682041650380972580064624519951808101508578950999055919377645323681448562195652215637418375701652544201727291448931246460221732662597172600081711");
    edwards97_Fq6::non_residue = edwards97_Fq("5");
    edwards97_Fq6::nqr = edwards97_Fq6(edwards97_Fq3(edwards97_Fq("3"),edwards97_Fq("0"),edwards97_Fq("0")),edwards97_Fq3::one());
    edwards97_Fq6::nqr_to_t = edwards97_Fq6(edwards97_Fq3::zero(),edwards97_Fq3(edwards97_Fq("0"),edwards97_Fq("242979017971813577150230832196"),edwards97_Fq("0")));
    edwards97_Fq6::Frobenius_coeffs_c1[0] = edwards97_Fq("1");
    edwards97_Fq6::Frobenius_coeffs_c1[1] = edwards97_Fq("394610617619779289602567691435");
    edwards97_Fq6::Frobenius_coeffs_c1[2] = edwards97_Fq("394610617619779289602567691434");
    edwards97_Fq6::Frobenius_coeffs_c1[3] = edwards97_Fq("565823376898968604518330826752");
    edwards97_Fq6::Frobenius_coeffs_c1[4] = edwards97_Fq("171212759279189314915763135318");
    edwards97_Fq6::Frobenius_coeffs_c1[5] = edwards97_Fq("171212759279189314915763135319");
    edwards97_Fq6::my_Fp2::non_residue = edwards97_Fq3::non_residue;
}

} // namespace libff
