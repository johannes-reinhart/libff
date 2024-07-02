#include <libff/algebra/curves/jubjub_ed97/jubjub_ed97_fields.hpp>

namespace libff {

bigint<jubjub_ed97_r_limbs> jubjub_ed97_modulus_r;
bigint<jubjub_ed97_q_limbs> jubjub_ed97_modulus_q;

void init_jubjub_ed97_fields()
{
    using bigint_r = bigint<jubjub_ed97_r_limbs>;
    using bigint_q = bigint<jubjub_ed97_q_limbs>;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    jubjub_ed97_modulus_r = bigint_r("17681980528092752387406989869");
    assert(jubjub_ed97_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        jubjub_ed97_Fr::Rsquared = bigint_r("9773283138541348269353355820");
        jubjub_ed97_Fr::Rcubed = bigint_r("9092061651246399697010690910");
        jubjub_ed97_Fr::inv = 0xb4f48634cb27aa5b;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        jubjub_ed97_Fr::Rsquared = bigint_r("7677772318542104462582833931");
        jubjub_ed97_Fr::Rcubed = bigint_r("4856437576305133042961554867");
        jubjub_ed97_Fr::inv = 0xcb27aa5b;
    }
    jubjub_ed97_Fr::num_bits = 94;
    jubjub_ed97_Fr::euler = bigint_r("8840990264046376193703494934");
    jubjub_ed97_Fr::s = 2;
    jubjub_ed97_Fr::t = bigint_r("4420495132023188096851747467");
    jubjub_ed97_Fr::t_minus_1_over_2 = bigint_r("2210247566011594048425873733");
    jubjub_ed97_Fr::multiplicative_generator = jubjub_ed97_Fr("6");
    jubjub_ed97_Fr::root_of_unity = jubjub_ed97_Fr("14545940802042798162886534517");
    jubjub_ed97_Fr::nqr = jubjub_ed97_Fr("2");
    jubjub_ed97_Fr::nqr_to_t = jubjub_ed97_Fr("3136039726049954224520455352");

    /* parameters for base field Fq */

    jubjub_ed97_modulus_q = bigint_q("141455844224742490147094691841");
    assert(jubjub_ed97_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        jubjub_ed97_Fq::Rsquared = bigint_q("75217608311129884939407662890");
        jubjub_ed97_Fq::Rcubed = bigint_q("128645697791147296497719755611");
        jubjub_ed97_Fq::inv = 0x69fd1db9f8ce7fff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        jubjub_ed97_Fq::Rsquared = bigint_q("75217608311129884939407662890");
        jubjub_ed97_Fq::Rcubed = bigint_q("128645697791147296497719755611");
        jubjub_ed97_Fq::inv = 0xf8ce7fff;
    }
    jubjub_ed97_Fq::num_bits = 97;
    jubjub_ed97_Fq::euler = bigint_q("70727922112371245073547345920");
    jubjub_ed97_Fq::s = 15;
    jubjub_ed97_Fq::t = bigint_q("4316889777366409001071005");
    jubjub_ed97_Fq::t_minus_1_over_2 = bigint_q("2158444888683204500535502");
    jubjub_ed97_Fq::multiplicative_generator = jubjub_ed97_Fq("13");
    jubjub_ed97_Fq::root_of_unity = jubjub_ed97_Fq("111466351910010104363354434347");
    jubjub_ed97_Fq::nqr = jubjub_ed97_Fq("13");
    jubjub_ed97_Fq::nqr_to_t = jubjub_ed97_Fq("111466351910010104363354434347");

    /* parameters for twist field Fq3 */

    jubjub_ed97_Fq3::euler = bigint<3*jubjub_ed97_q_limbs>("1415248454330950213824664324088701009430335730854515614315458818484823798747890296668160");
    jubjub_ed97_Fq3::s = 15;
    jubjub_ed97_Fq3::t = bigint<3*jubjub_ed97_q_limbs>("86379910542660535511759297124554504970113264822663306537808765776661608810296038615");
    jubjub_ed97_Fq3::t_minus_1_over_2 = bigint<3*jubjub_ed97_q_limbs>("43189955271330267755879648562277252485056632411331653268904382888330804405148019307");
    jubjub_ed97_Fq3::non_residue = jubjub_ed97_Fq("13");
    jubjub_ed97_Fq3::nqr = jubjub_ed97_Fq3(jubjub_ed97_Fq("13"),jubjub_ed97_Fq("0"),jubjub_ed97_Fq("0"));
    jubjub_ed97_Fq3::nqr_to_t = jubjub_ed97_Fq3(jubjub_ed97_Fq("41627948615973198239244211885"),jubjub_ed97_Fq("0"),jubjub_ed97_Fq("0"));
    jubjub_ed97_Fq3::Frobenius_coeffs_c1[0] = jubjub_ed97_Fq("1");
    jubjub_ed97_Fq3::Frobenius_coeffs_c1[1] = jubjub_ed97_Fq("141455844224741134077046751229");
    jubjub_ed97_Fq3::Frobenius_coeffs_c1[2] = jubjub_ed97_Fq("1356070047940611");
    jubjub_ed97_Fq3::Frobenius_coeffs_c2[0] = jubjub_ed97_Fq("1");
    jubjub_ed97_Fq3::Frobenius_coeffs_c2[1] = jubjub_ed97_Fq("1356070047940611");
    jubjub_ed97_Fq3::Frobenius_coeffs_c2[2] = jubjub_ed97_Fq("141455844224741134077046751229");

    /* parameters for Fq6 */

    jubjub_ed97_Fq6::euler = bigint<6*jubjub_ed97_q_limbs>("4005856374972287346084358535003138324060857269877739947402181403384023057462738434804023761723442065397649857465418089797361316782026939219286663847853068835094249374908907520");
    jubjub_ed97_Fq6::s = 16;
    jubjub_ed97_Fq6::t = bigint<6*jubjub_ed97_q_limbs>("122249034880746073794078324432468820924708778987968138043279461773194063032920484460572014212751527874684138716596011041179239403748380713479207270747469141695991497037015");
    jubjub_ed97_Fq6::t_minus_1_over_2 = bigint<6*jubjub_ed97_q_limbs>("61124517440373036897039162216234410462354389493984069021639730886597031516460242230286007106375763937342069358298005520589619701874190356739603635373734570847995748518507");
    jubjub_ed97_Fq6::non_residue = jubjub_ed97_Fq("13");
    jubjub_ed97_Fq6::nqr = jubjub_ed97_Fq6(jubjub_ed97_Fq3(jubjub_ed97_Fq("4"),jubjub_ed97_Fq("0"),jubjub_ed97_Fq("0")),jubjub_ed97_Fq3::one());
    jubjub_ed97_Fq6::nqr_to_t = jubjub_ed97_Fq6(jubjub_ed97_Fq3::zero(),jubjub_ed97_Fq3(jubjub_ed97_Fq("0"),jubjub_ed97_Fq("131422745497037914015619005162"),jubjub_ed97_Fq("0")));
    jubjub_ed97_Fq6::Frobenius_coeffs_c1[0] = jubjub_ed97_Fq("1");
    jubjub_ed97_Fq6::Frobenius_coeffs_c1[1] = jubjub_ed97_Fq("141455844224741134077046751230");
    jubjub_ed97_Fq6::Frobenius_coeffs_c1[2] = jubjub_ed97_Fq("141455844224741134077046751229");
    jubjub_ed97_Fq6::Frobenius_coeffs_c1[3] = jubjub_ed97_Fq("141455844224742490147094691840");
    jubjub_ed97_Fq6::Frobenius_coeffs_c1[4] = jubjub_ed97_Fq("1356070047940611");
    jubjub_ed97_Fq6::Frobenius_coeffs_c1[5] = jubjub_ed97_Fq("1356070047940612");
    jubjub_ed97_Fq6::my_Fp2::non_residue = jubjub_ed97_Fq3::non_residue;
}

} // namespace libff
