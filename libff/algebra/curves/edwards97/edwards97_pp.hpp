#ifndef EDWARDS97_PP_HPP_
#define EDWARDS97_PP_HPP_
#include <libff/algebra/curves/edwards97/edwards97_g1.hpp>
#include <libff/algebra/curves/edwards97/edwards97_g2.hpp>
#include <libff/algebra/curves/edwards97/edwards97_init.hpp>
#include <libff/algebra/curves/edwards97/edwards97_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class edwards97_pp {
public:
    typedef edwards97_Fr Fp_type;
    typedef edwards97_G1 G1_type;
    typedef edwards97_G2 G2_type;
    typedef edwards97_G1_precomp G1_precomp_type;
    typedef edwards97_G2_precomp G2_precomp_type;
    typedef edwards97_Fq Fq_type;
    typedef edwards97_Fq3 Fqe_type;
    typedef edwards97_Fq6 Fqk_type;
    typedef edwards97_GT GT_type;

    static const bool has_affine_pairing = false;

    static void init_public_params();
    static edwards97_GT final_exponentiation(const edwards97_Fq6 &elt);
    static edwards97_G1_precomp precompute_G1(const edwards97_G1 &P);
    static edwards97_G2_precomp precompute_G2(const edwards97_G2 &Q);
    static edwards97_Fq6 miller_loop(const edwards97_G1_precomp &prec_P,
                                   const edwards97_G2_precomp &prec_Q);
    static edwards97_Fq6 double_miller_loop(const edwards97_G1_precomp &prec_P1,
                                          const edwards97_G2_precomp &prec_Q1,
                                          const edwards97_G1_precomp &prec_P2,
                                          const edwards97_G2_precomp &prec_Q2);
    /* the following are used in test files */
    static edwards97_Fq6 pairing(const edwards97_G1 &P,
                               const edwards97_G2 &Q);
    static edwards97_Fq6 reduced_pairing(const edwards97_G1 &P,
                                       const edwards97_G2 &Q);
};

} // namespace libff
#endif // EDWARDS97_PP_HPP_
