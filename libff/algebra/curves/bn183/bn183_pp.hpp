#ifndef BN183_PP_HPP_
#define BN183_PP_HPP_
#include <libff/algebra/curves/bn183/bn183_g1.hpp>
#include <libff/algebra/curves/bn183/bn183_g2.hpp>
#include <libff/algebra/curves/bn183/bn183_init.hpp>
#include <libff/algebra/curves/bn183/bn183_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class bn183_pp {
public:
    typedef bn183_Fr Fp_type;
    typedef bn183_G1 G1_type;
    typedef bn183_G2 G2_type;
    typedef bn183_G1_precomp G1_precomp_type;
    typedef bn183_G2_precomp G2_precomp_type;
    typedef bn183_Fq Fq_type;
    typedef bn183_Fq2 Fqe_type;
    typedef bn183_Fq12 Fqk_type;
    typedef bn183_GT GT_type;

    static const bool has_affine_pairing = false;

    static void init_public_params();
    static bn183_GT final_exponentiation(const bn183_Fq12 &elt);
    static bn183_G1_precomp precompute_G1(const bn183_G1 &P);
    static bn183_G2_precomp precompute_G2(const bn183_G2 &Q);
    static bn183_Fq12 miller_loop(const bn183_G1_precomp &prec_P,
                                      const bn183_G2_precomp &prec_Q);
    static bn183_Fq12 double_miller_loop(const bn183_G1_precomp &prec_P1,
                                             const bn183_G2_precomp &prec_Q1,
                                             const bn183_G1_precomp &prec_P2,
                                             const bn183_G2_precomp &prec_Q2);
    static bn183_Fq12 pairing(const bn183_G1 &P,
                                  const bn183_G2 &Q);
    static bn183_Fq12 reduced_pairing(const bn183_G1 &P,
                                          const bn183_G2 &Q);
};

} // namespace libff

#endif // BN183_PP_HPP_
