#include <libff/algebra/curves/edwards97/edwards97_pp.hpp>

namespace libff {

void edwards97_pp::init_public_params()
{
    init_edwards97_params();
}

edwards97_GT edwards97_pp::final_exponentiation(const edwards97_Fq6 &elt)
{
    return edwards97_final_exponentiation(elt);
}

edwards97_G1_precomp edwards97_pp::precompute_G1(const edwards97_G1 &P)
{
    return edwards97_precompute_G1(P);
}

edwards97_G2_precomp edwards97_pp::precompute_G2(const edwards97_G2 &Q)
{
    return edwards97_precompute_G2(Q);
}

edwards97_Fq6 edwards97_pp::miller_loop(const edwards97_G1_precomp &prec_P,
                                    const edwards97_G2_precomp &prec_Q)
{
    return edwards97_miller_loop(prec_P, prec_Q);
}

edwards97_Fq6 edwards97_pp::double_miller_loop(const edwards97_G1_precomp &prec_P1,
                                           const edwards97_G2_precomp &prec_Q1,
                                           const edwards97_G1_precomp &prec_P2,
                                           const edwards97_G2_precomp &prec_Q2)
{
    return edwards97_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

edwards97_Fq6 edwards97_pp::pairing(const edwards97_G1 &P,
                                const edwards97_G2 &Q)
{
    return edwards97_pairing(P, Q);
}

edwards97_Fq6 edwards97_pp::reduced_pairing(const edwards97_G1 &P,
                                        const edwards97_G2 &Q)
{
    return edwards97_reduced_pairing(P, Q);
}

} // namespace libff
