#include <libff/algebra/curves/edwards58/edwards58_fields.hpp>

namespace libff {

bigint<edwards58_r_limbs> edwards58_modulus_r;
bigint<edwards58_q_limbs> edwards58_modulus_q;

void init_edwards58_fields()
{
    using bigint_r = bigint<edwards58_r_limbs>;
    using bigint_q = bigint<edwards58_q_limbs>;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    edwards58_modulus_r = bigint_r("211006452744585217");
    assert(edwards58_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        edwards58_Fr::Rsquared = bigint_r("69518267247397641");
        edwards58_Fr::Rcubed = bigint_r("193528881679320662");
        edwards58_Fr::inv = 0x8de40f003527ffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        edwards58_Fr::Rsquared = bigint_r("69518267247397641");
        edwards58_Fr::Rcubed = bigint_r("193528881679320662");
        edwards58_Fr::inv = 0x3527ffff;
    }
    edwards58_Fr::num_bits = 58;
    edwards58_Fr::euler = bigint_r("105503226372292608");
    edwards58_Fr::s = 19;
    edwards58_Fr::t = bigint_r("402462869157");
    edwards58_Fr::t_minus_1_over_2 = bigint_r("201231434578");
    edwards58_Fr::multiplicative_generator = edwards58_Fr("7");
    edwards58_Fr::root_of_unity = edwards58_Fr("101722017355283520");
    edwards58_Fr::nqr = edwards58_Fr("5");
    edwards58_Fr::nqr_to_t = edwards58_Fr("85346003589280634");

    /* parameters for base field Fq */

    edwards58_modulus_q = bigint_q("844025809322115073");
    assert(edwards58_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        edwards58_Fq::Rsquared = bigint_q("81106780184426238");
        edwards58_Fq::Rcubed = bigint_q("231719045362174245");
        edwards58_Fq::inv = 0x1307f2c071e7ffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        edwards58_Fq::Rsquared = bigint_q("81106780184426238");
        edwards58_Fq::Rcubed = bigint_q("231719045362174245");
        edwards58_Fq::inv = 0x71e7ffff;
    }
    edwards58_Fq::num_bits = 60;
    edwards58_Fq::euler = bigint_q("422012904661057536");
    edwards58_Fq::s = 19;
    edwards58_Fq::t = bigint_q("1609851473469");
    edwards58_Fq::t_minus_1_over_2 = bigint_q("804925736734");
    edwards58_Fq::multiplicative_generator = edwards58_Fq("10");
    edwards58_Fq::root_of_unity = edwards58_Fq("749735321661346398");
    edwards58_Fq::nqr = edwards58_Fq("5");
    edwards58_Fq::nqr_to_t = edwards58_Fq("473066918512820194");

    /* parameters for twist field Fq3 */

    edwards58_Fq3::euler = bigint<3*edwards58_q_limbs>("300633370207235162806043354206417528559416271963947008");
    edwards58_Fq3::s = 19;
    edwards58_Fq3::t = bigint<3*edwards58_q_limbs>("1146825295285168315147565285516424288022675598007");
    edwards58_Fq3::t_minus_1_over_2 = bigint<3*edwards58_q_limbs>("573412647642584157573782642758212144011337799003");
    edwards58_Fq3::non_residue = edwards58_Fq("10");
    edwards58_Fq3::nqr = edwards58_Fq3(edwards58_Fq("5"),edwards58_Fq("0"),edwards58_Fq("0"));
    edwards58_Fq3::nqr_to_t = edwards58_Fq3(edwards58_Fq("628058373490290683"),edwards58_Fq("0"),edwards58_Fq("0"));
    edwards58_Fq3::Frobenius_coeffs_c1[0] = edwards58_Fq("1");
    edwards58_Fq3::Frobenius_coeffs_c1[1] = edwards58_Fq("126784148197437870");
    edwards58_Fq3::Frobenius_coeffs_c1[2] = edwards58_Fq("717241661124677202");
    edwards58_Fq3::Frobenius_coeffs_c2[0] = edwards58_Fq("1");
    edwards58_Fq3::Frobenius_coeffs_c2[1] = edwards58_Fq("717241661124677202");
    edwards58_Fq3::Frobenius_coeffs_c2[2] = edwards58_Fq("126784148197437870");

    /* parameters for Fq6 */

    edwards58_Fq6::euler = bigint<6*edwards58_q_limbs>("180760846564321021593410705171210498602927285896697201913540158600790427791413753146341399605088884392198144");
    edwards58_Fq6::s = 20;
    edwards58_Fq6::t = bigint<6*edwards58_q_limbs>("344773953560487788378545198767109868245939800065416721179085080339032035429790025990183638773134011063");
    edwards58_Fq6::t_minus_1_over_2 = bigint<6*edwards58_q_limbs>("172386976780243894189272599383554934122969900032708360589542540169516017714895012995091819386567005531");
    edwards58_Fq6::non_residue = edwards58_Fq("10");
    edwards58_Fq6::nqr = edwards58_Fq6(edwards58_Fq3(edwards58_Fq("3"),edwards58_Fq("0"),edwards58_Fq("0")),edwards58_Fq3::one());
    edwards58_Fq6::nqr_to_t = edwards58_Fq6(edwards58_Fq3::zero(),edwards58_Fq3(edwards58_Fq("0"),edwards58_Fq("587312879757227502"),edwards58_Fq("0")));
    edwards58_Fq6::Frobenius_coeffs_c1[0] = edwards58_Fq("1");
    edwards58_Fq6::Frobenius_coeffs_c1[1] = edwards58_Fq("126784148197437871");
    edwards58_Fq6::Frobenius_coeffs_c1[2] = edwards58_Fq("126784148197437870");
    edwards58_Fq6::Frobenius_coeffs_c1[3] = edwards58_Fq("844025809322115072");
    edwards58_Fq6::Frobenius_coeffs_c1[4] = edwards58_Fq("717241661124677202");
    edwards58_Fq6::Frobenius_coeffs_c1[5] = edwards58_Fq("717241661124677203");
    edwards58_Fq6::my_Fp2::non_residue = edwards58_Fq3::non_residue;
}

} // namespace libff
