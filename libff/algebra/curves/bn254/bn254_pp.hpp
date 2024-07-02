#ifndef BN254_PP_HPP_
#define BN254_PP_HPP_
#include <libff/algebra/curves/bn254/bn254_g1.hpp>
#include <libff/algebra/curves/bn254/bn254_g2.hpp>
#include <libff/algebra/curves/bn254/bn254_init.hpp>
#include <libff/algebra/curves/bn254/bn254_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class bn254_pp {
public:
    typedef bn254_Fr Fp_type;
    typedef bn254_G1 G1_type;
    typedef bn254_G2 G2_type;
    typedef bn254_G1_precomp G1_precomp_type;
    typedef bn254_G2_precomp G2_precomp_type;
    typedef bn254_Fq Fq_type;
    typedef bn254_Fq2 Fqe_type;
    typedef bn254_Fq12 Fqk_type;
    typedef bn254_GT GT_type;

    static const bool has_affine_pairing = false;

    static void init_public_params();
    static bn254_GT final_exponentiation(const bn254_Fq12 &elt);
    static bn254_G1_precomp precompute_G1(const bn254_G1 &P);
    static bn254_G2_precomp precompute_G2(const bn254_G2 &Q);
    static bn254_Fq12 miller_loop(const bn254_G1_precomp &prec_P,
                                      const bn254_G2_precomp &prec_Q);
    static bn254_Fq12 double_miller_loop(const bn254_G1_precomp &prec_P1,
                                             const bn254_G2_precomp &prec_Q1,
                                             const bn254_G1_precomp &prec_P2,
                                             const bn254_G2_precomp &prec_Q2);
    static bn254_Fq12 pairing(const bn254_G1 &P,
                                  const bn254_G2 &Q);
    static bn254_Fq12 reduced_pairing(const bn254_G1 &P,
                                          const bn254_G2 &Q);
};

} // namespace libff

#endif // BN254_PP_HPP_
