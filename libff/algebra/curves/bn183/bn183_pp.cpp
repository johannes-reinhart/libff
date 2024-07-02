#include <libff/algebra/curves/bn183/bn183_pp.hpp>

namespace libff {

void bn183_pp::init_public_params()
{
    init_bn183_params();
}

bn183_GT bn183_pp::final_exponentiation(const bn183_Fq12 &elt)
{
    return bn183_final_exponentiation(elt);
}

bn183_G1_precomp bn183_pp::precompute_G1(const bn183_G1 &P)
{
    return bn183_precompute_G1(P);
}

bn183_G2_precomp bn183_pp::precompute_G2(const bn183_G2 &Q)
{
    return bn183_precompute_G2(Q);
}

bn183_Fq12 bn183_pp::miller_loop(const bn183_G1_precomp &prec_P,
                                         const bn183_G2_precomp &prec_Q)
{
    return bn183_miller_loop(prec_P, prec_Q);
}

bn183_Fq12 bn183_pp::double_miller_loop(const bn183_G1_precomp &prec_P1,
                                                const bn183_G2_precomp &prec_Q1,
                                                const bn183_G1_precomp &prec_P2,
                                                const bn183_G2_precomp &prec_Q2)
{
    return bn183_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

bn183_Fq12 bn183_pp::pairing(const bn183_G1 &P,
                                     const bn183_G2 &Q)
{
    return bn183_pairing(P, Q);
}

bn183_Fq12 bn183_pp::reduced_pairing(const bn183_G1 &P,
                                             const bn183_G2 &Q)
{
    return bn183_reduced_pairing(P, Q);
}

} // namespace libff
