#include <libff/algebra/curves/edwards61/edwards61_fields.hpp>

namespace libff {

bigint<edwards61_r_limbs> edwards61_modulus_r;
bigint<edwards61_q_limbs> edwards61_modulus_q;

void init_edwards61_fields()
{
    using bigint_r = bigint<edwards61_r_limbs>;
    using bigint_q = bigint<edwards61_q_limbs>;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    edwards61_modulus_r = bigint_r("1235486033917771777");
    assert(edwards61_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        edwards61_Fr::Rsquared = bigint_r("1073255652645550202");
        edwards61_Fr::Rcubed = bigint_r("357882241080121889");
        edwards61_Fr::inv = 0x5084f000809fffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        edwards61_Fr::Rsquared = bigint_r("1073255652645550202");
        edwards61_Fr::Rcubed = bigint_r("357882241080121889");
        edwards61_Fr::inv = 0x809fffff;
    }
    edwards61_Fr::num_bits = 61;
    edwards61_Fr::euler = bigint_r("617743016958885888");
    edwards61_Fr::s = 21;
    edwards61_Fr::t = bigint_r("589125649413");
    edwards61_Fr::t_minus_1_over_2 = bigint_r("294562824706");
    edwards61_Fr::multiplicative_generator = edwards61_Fr("13");
    edwards61_Fr::root_of_unity = edwards61_Fr("909855932006627559");
    edwards61_Fr::nqr = edwards61_Fr("5");
    edwards61_Fr::nqr_to_t = edwards61_Fr("468172926162659554");

    /* parameters for base field Fq */

    edwards61_modulus_q = bigint_q("4941944131663429633");
    assert(edwards61_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        edwards61_Fq::Rsquared = bigint_q("1307098666210323107");
        edwards61_Fq::Rcubed = bigint_q("620669261281184929");
        edwards61_Fq::inv = 0x9bd42c01139fffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        edwards61_Fq::Rsquared = bigint_q("1307098666210323107");
        edwards61_Fq::Rcubed = bigint_q("620669261281184929");
        edwards61_Fq::inv = 0x139fffff;
    }
    edwards61_Fq::num_bits = 63;
    edwards61_Fq::euler = bigint_q("2470972065831714816");
    edwards61_Fq::s = 21;
    edwards61_Fq::t = bigint_q("2356502595741");
    edwards61_Fq::t_minus_1_over_2 = bigint_q("1178251297870");
    edwards61_Fq::multiplicative_generator = edwards61_Fq("5");
    edwards61_Fq::root_of_unity = edwards61_Fq("4812765192829825337");
    edwards61_Fq::nqr = edwards61_Fq("5");
    edwards61_Fq::nqr_to_t = edwards61_Fq("4812765192829825337");

    /* parameters for twist field Fq3 */

    edwards61_Fq3::euler = bigint<3*edwards61_q_limbs>("60348085728057696086447257716996238190658477546110189568");
    edwards61_Fq3::s = 21;
    edwards61_Fq3::t = bigint<3*edwards61_q_limbs>("57552419403131195150801904408451307478579022928343");
    edwards61_Fq3::t_minus_1_over_2 = bigint<3*edwards61_q_limbs>("28776209701565597575400952204225653739289511464171");
    edwards61_Fq3::non_residue = edwards61_Fq("5");
    edwards61_Fq3::nqr = edwards61_Fq3(edwards61_Fq("5"),edwards61_Fq("0"),edwards61_Fq("0"));
    edwards61_Fq3::nqr_to_t = edwards61_Fq3(edwards61_Fq("4129413696311215253"),edwards61_Fq("0"),edwards61_Fq("0"));
    edwards61_Fq3::Frobenius_coeffs_c1[0] = edwards61_Fq("1");
    edwards61_Fq3::Frobenius_coeffs_c1[1] = edwards61_Fq("875438387576607505");
    edwards61_Fq3::Frobenius_coeffs_c1[2] = edwards61_Fq("4066505744086822127");
    edwards61_Fq3::Frobenius_coeffs_c2[0] = edwards61_Fq("1");
    edwards61_Fq3::Frobenius_coeffs_c2[1] = edwards61_Fq("4066505744086822127");
    edwards61_Fq3::Frobenius_coeffs_c2[2] = edwards61_Fq("875438387576607505");

    /* parameters for Fq6 */

    edwards61_Fq6::euler = bigint<6*edwards61_q_limbs>("7283782902082001973452363488299976838804677640302607409910163597201576981962574354726579877718666999086012432384");
    edwards61_Fq6::s = 22;
    edwards61_Fq6::t = bigint<6*edwards61_q_limbs>("3473178339997292505956823105001438540842379398490241723017770575142658701878821542132654131755193233054167");
    edwards61_Fq6::t_minus_1_over_2 = bigint<6*edwards61_q_limbs>("1736589169998646252978411552500719270421189699245120861508885287571329350939410771066327065877596616527083");
    edwards61_Fq6::non_residue = edwards61_Fq("5");
    edwards61_Fq6::nqr = edwards61_Fq6(edwards61_Fq3(edwards61_Fq("2"),edwards61_Fq("0"),edwards61_Fq("0")),edwards61_Fq3::one());
    edwards61_Fq6::nqr_to_t = edwards61_Fq6(edwards61_Fq3::zero(),edwards61_Fq3(edwards61_Fq("0"),edwards61_Fq("3520480479100147022"),edwards61_Fq("0")));
    edwards61_Fq6::Frobenius_coeffs_c1[0] = edwards61_Fq("1");
    edwards61_Fq6::Frobenius_coeffs_c1[1] = edwards61_Fq("875438387576607506");
    edwards61_Fq6::Frobenius_coeffs_c1[2] = edwards61_Fq("875438387576607505");
    edwards61_Fq6::Frobenius_coeffs_c1[3] = edwards61_Fq("4941944131663429632");
    edwards61_Fq6::Frobenius_coeffs_c1[4] = edwards61_Fq("4066505744086822127");
    edwards61_Fq6::Frobenius_coeffs_c1[5] = edwards61_Fq("4066505744086822128");
    edwards61_Fq6::my_Fp2::non_residue = edwards61_Fq3::non_residue;
}

} // namespace libff
