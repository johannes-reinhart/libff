#ifndef BN124_PP_HPP_
#define BN124_PP_HPP_
#include <libff/algebra/curves/bn124/bn124_g1.hpp>
#include <libff/algebra/curves/bn124/bn124_g2.hpp>
#include <libff/algebra/curves/bn124/bn124_init.hpp>
#include <libff/algebra/curves/bn124/bn124_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class bn124_pp {
public:
    typedef bn124_Fr Fp_type;
    typedef bn124_G1 G1_type;
    typedef bn124_G2 G2_type;
    typedef bn124_G1_precomp G1_precomp_type;
    typedef bn124_G2_precomp G2_precomp_type;
    typedef bn124_Fq Fq_type;
    typedef bn124_Fq2 Fqe_type;
    typedef bn124_Fq12 Fqk_type;
    typedef bn124_GT GT_type;

    static const bool has_affine_pairing = false;

    static void init_public_params();
    static bn124_GT final_exponentiation(const bn124_Fq12 &elt);
    static bn124_G1_precomp precompute_G1(const bn124_G1 &P);
    static bn124_G2_precomp precompute_G2(const bn124_G2 &Q);
    static bn124_Fq12 miller_loop(const bn124_G1_precomp &prec_P,
                                      const bn124_G2_precomp &prec_Q);
    static bn124_Fq12 double_miller_loop(const bn124_G1_precomp &prec_P1,
                                             const bn124_G2_precomp &prec_Q1,
                                             const bn124_G1_precomp &prec_P2,
                                             const bn124_G2_precomp &prec_Q2);
    static bn124_Fq12 pairing(const bn124_G1 &P,
                                  const bn124_G2 &Q);
    static bn124_Fq12 reduced_pairing(const bn124_G1 &P,
                                          const bn124_G2 &Q);
};

} // namespace libff

#endif // BN124_PP_HPP_
