#include <libff/algebra/curves/bn124/bn124_pp.hpp>

namespace libff {

void bn124_pp::init_public_params()
{
    init_bn124_params();
}

bn124_GT bn124_pp::final_exponentiation(const bn124_Fq12 &elt)
{
    return bn124_final_exponentiation(elt);
}

bn124_G1_precomp bn124_pp::precompute_G1(const bn124_G1 &P)
{
    return bn124_precompute_G1(P);
}

bn124_G2_precomp bn124_pp::precompute_G2(const bn124_G2 &Q)
{
    return bn124_precompute_G2(Q);
}

bn124_Fq12 bn124_pp::miller_loop(const bn124_G1_precomp &prec_P,
                                         const bn124_G2_precomp &prec_Q)
{
    return bn124_miller_loop(prec_P, prec_Q);
}

bn124_Fq12 bn124_pp::double_miller_loop(const bn124_G1_precomp &prec_P1,
                                                const bn124_G2_precomp &prec_Q1,
                                                const bn124_G1_precomp &prec_P2,
                                                const bn124_G2_precomp &prec_Q2)
{
    return bn124_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

bn124_Fq12 bn124_pp::pairing(const bn124_G1 &P,
                                     const bn124_G2 &Q)
{
    return bn124_pairing(P, Q);
}

bn124_Fq12 bn124_pp::reduced_pairing(const bn124_G1 &P,
                                             const bn124_G2 &Q)
{
    return bn124_reduced_pairing(P, Q);
}

} // namespace libff
