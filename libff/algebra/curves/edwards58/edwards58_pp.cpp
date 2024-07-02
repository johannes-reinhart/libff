#include <libff/algebra/curves/edwards58/edwards58_pp.hpp>

namespace libff {

void edwards58_pp::init_public_params()
{
    init_edwards58_params();
}

edwards58_GT edwards58_pp::final_exponentiation(const edwards58_Fq6 &elt)
{
    return edwards58_final_exponentiation(elt);
}

edwards58_G1_precomp edwards58_pp::precompute_G1(const edwards58_G1 &P)
{
    return edwards58_precompute_G1(P);
}

edwards58_G2_precomp edwards58_pp::precompute_G2(const edwards58_G2 &Q)
{
    return edwards58_precompute_G2(Q);
}

edwards58_Fq6 edwards58_pp::miller_loop(const edwards58_G1_precomp &prec_P,
                                    const edwards58_G2_precomp &prec_Q)
{
    return edwards58_miller_loop(prec_P, prec_Q);
}

edwards58_Fq6 edwards58_pp::double_miller_loop(const edwards58_G1_precomp &prec_P1,
                                           const edwards58_G2_precomp &prec_Q1,
                                           const edwards58_G1_precomp &prec_P2,
                                           const edwards58_G2_precomp &prec_Q2)
{
    return edwards58_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

edwards58_Fq6 edwards58_pp::pairing(const edwards58_G1 &P,
                                const edwards58_G2 &Q)
{
    return edwards58_pairing(P, Q);
}

edwards58_Fq6 edwards58_pp::reduced_pairing(const edwards58_G1 &P,
                                        const edwards58_G2 &Q)
{
    return edwards58_reduced_pairing(P, Q);
}

} // namespace libff
