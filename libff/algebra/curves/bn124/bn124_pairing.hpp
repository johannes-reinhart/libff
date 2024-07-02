#ifndef BN124_PAIRING_HPP_
#define BN124_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/bn124/bn124_init.hpp>

namespace libff {

/* final exponentiation */

bn124_GT bn124_final_exponentiation(const bn124_Fq12 &elt);

/* ate pairing */

struct bn124_ate_G1_precomp {
    bn124_Fq PX;
    bn124_Fq PY;

    bool operator==(const bn124_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn124_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, bn124_ate_G1_precomp &prec_P);
};

struct bn124_ate_ell_coeffs {
    bn124_Fq2 ell_0;
    bn124_Fq2 ell_VW;
    bn124_Fq2 ell_VV;

    bool operator==(const bn124_ate_ell_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn124_ate_ell_coeffs &c);
    friend std::istream& operator>>(std::istream &in, bn124_ate_ell_coeffs &c);
};

struct bn124_ate_G2_precomp {
    bn124_Fq2 QX;
    bn124_Fq2 QY;
    std::vector<bn124_ate_ell_coeffs> coeffs;

    bool operator==(const bn124_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn124_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, bn124_ate_G2_precomp &prec_Q);
};

bn124_ate_G1_precomp bn124_ate_precompute_G1(const bn124_G1& P);
bn124_ate_G2_precomp bn124_ate_precompute_G2(const bn124_G2& Q);

bn124_Fq12 bn124_ate_miller_loop(const bn124_ate_G1_precomp &prec_P,
                              const bn124_ate_G2_precomp &prec_Q);
bn124_Fq12 bn124_ate_double_miller_loop(const bn124_ate_G1_precomp &prec_P1,
                                     const bn124_ate_G2_precomp &prec_Q1,
                                     const bn124_ate_G1_precomp &prec_P2,
                                     const bn124_ate_G2_precomp &prec_Q2);

bn124_Fq12 bn124_ate_pairing(const bn124_G1& P,
                          const bn124_G2 &Q);
bn124_GT bn124_ate_reduced_pairing(const bn124_G1 &P,
                                 const bn124_G2 &Q);

/* choice of pairing */

typedef bn124_ate_G1_precomp bn124_G1_precomp;
typedef bn124_ate_G2_precomp bn124_G2_precomp;

bn124_G1_precomp bn124_precompute_G1(const bn124_G1& P);

bn124_G2_precomp bn124_precompute_G2(const bn124_G2& Q);

bn124_Fq12 bn124_miller_loop(const bn124_G1_precomp &prec_P,
                          const bn124_G2_precomp &prec_Q);

bn124_Fq12 bn124_double_miller_loop(const bn124_G1_precomp &prec_P1,
                                 const bn124_G2_precomp &prec_Q1,
                                 const bn124_G1_precomp &prec_P2,
                                 const bn124_G2_precomp &prec_Q2);

bn124_Fq12 bn124_pairing(const bn124_G1& P,
                      const bn124_G2 &Q);

bn124_GT bn124_reduced_pairing(const bn124_G1 &P,
                             const bn124_G2 &Q);

bn124_GT bn124_affine_reduced_pairing(const bn124_G1 &P,
                                    const bn124_G2 &Q);

} // namespace libff
#endif // BN124_PAIRING_HPP_
