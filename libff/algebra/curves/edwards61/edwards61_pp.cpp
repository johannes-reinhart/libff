#include <libff/algebra/curves/edwards61/edwards61_pp.hpp>

namespace libff {

void edwards61_pp::init_public_params()
{
    init_edwards61_params();
}

edwards61_GT edwards61_pp::final_exponentiation(const edwards61_Fq6 &elt)
{
    return edwards61_final_exponentiation(elt);
}

edwards61_G1_precomp edwards61_pp::precompute_G1(const edwards61_G1 &P)
{
    return edwards61_precompute_G1(P);
}

edwards61_G2_precomp edwards61_pp::precompute_G2(const edwards61_G2 &Q)
{
    return edwards61_precompute_G2(Q);
}

edwards61_Fq6 edwards61_pp::miller_loop(const edwards61_G1_precomp &prec_P,
                                    const edwards61_G2_precomp &prec_Q)
{
    return edwards61_miller_loop(prec_P, prec_Q);
}

edwards61_Fq6 edwards61_pp::double_miller_loop(const edwards61_G1_precomp &prec_P1,
                                           const edwards61_G2_precomp &prec_Q1,
                                           const edwards61_G1_precomp &prec_P2,
                                           const edwards61_G2_precomp &prec_Q2)
{
    return edwards61_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

edwards61_Fq6 edwards61_pp::pairing(const edwards61_G1 &P,
                                const edwards61_G2 &Q)
{
    return edwards61_pairing(P, Q);
}

edwards61_Fq6 edwards61_pp::reduced_pairing(const edwards61_G1 &P,
                                        const edwards61_G2 &Q)
{
    return edwards61_reduced_pairing(P, Q);
}

} // namespace libff
