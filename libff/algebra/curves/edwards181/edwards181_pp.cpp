#include <libff/algebra/curves/edwards181/edwards181_pp.hpp>

namespace libff {

void edwards181_pp::init_public_params()
{
    init_edwards181_params();
}

edwards181_GT edwards181_pp::final_exponentiation(const edwards181_Fq6 &elt)
{
    return edwards181_final_exponentiation(elt);
}

edwards181_G1_precomp edwards181_pp::precompute_G1(const edwards181_G1 &P)
{
    return edwards181_precompute_G1(P);
}

edwards181_G2_precomp edwards181_pp::precompute_G2(const edwards181_G2 &Q)
{
    return edwards181_precompute_G2(Q);
}

edwards181_Fq6 edwards181_pp::miller_loop(const edwards181_G1_precomp &prec_P,
                                    const edwards181_G2_precomp &prec_Q)
{
    return edwards181_miller_loop(prec_P, prec_Q);
}

edwards181_Fq6 edwards181_pp::double_miller_loop(const edwards181_G1_precomp &prec_P1,
                                           const edwards181_G2_precomp &prec_Q1,
                                           const edwards181_G1_precomp &prec_P2,
                                           const edwards181_G2_precomp &prec_Q2)
{
    return edwards181_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

edwards181_Fq6 edwards181_pp::pairing(const edwards181_G1 &P,
                                const edwards181_G2 &Q)
{
    return edwards181_pairing(P, Q);
}

edwards181_Fq6 edwards181_pp::reduced_pairing(const edwards181_G1 &P,
                                        const edwards181_G2 &Q)
{
    return edwards181_reduced_pairing(P, Q);
}

} // namespace libff
