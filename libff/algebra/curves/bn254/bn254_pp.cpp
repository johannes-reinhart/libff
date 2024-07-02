#include <libff/algebra/curves/bn254/bn254_pp.hpp>

namespace libff {

void bn254_pp::init_public_params()
{
    init_bn254_params();
}

bn254_GT bn254_pp::final_exponentiation(const bn254_Fq12 &elt)
{
    return bn254_final_exponentiation(elt);
}

bn254_G1_precomp bn254_pp::precompute_G1(const bn254_G1 &P)
{
    return bn254_precompute_G1(P);
}

bn254_G2_precomp bn254_pp::precompute_G2(const bn254_G2 &Q)
{
    return bn254_precompute_G2(Q);
}

bn254_Fq12 bn254_pp::miller_loop(const bn254_G1_precomp &prec_P,
                                         const bn254_G2_precomp &prec_Q)
{
    return bn254_miller_loop(prec_P, prec_Q);
}

bn254_Fq12 bn254_pp::double_miller_loop(const bn254_G1_precomp &prec_P1,
                                                const bn254_G2_precomp &prec_Q1,
                                                const bn254_G1_precomp &prec_P2,
                                                const bn254_G2_precomp &prec_Q2)
{
    return bn254_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

bn254_Fq12 bn254_pp::pairing(const bn254_G1 &P,
                                     const bn254_G2 &Q)
{
    return bn254_pairing(P, Q);
}

bn254_Fq12 bn254_pp::reduced_pairing(const bn254_G1 &P,
                                             const bn254_G2 &Q)
{
    return bn254_reduced_pairing(P, Q);
}

} // namespace libff
