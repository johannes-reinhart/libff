#ifndef BN183_PAIRING_HPP_
#define BN183_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/bn183/bn183_init.hpp>

namespace libff {

/* final exponentiation */

bn183_GT bn183_final_exponentiation(const bn183_Fq12 &elt);

/* ate pairing */

struct bn183_ate_G1_precomp {
    bn183_Fq PX;
    bn183_Fq PY;

    bool operator==(const bn183_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn183_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, bn183_ate_G1_precomp &prec_P);
};

struct bn183_ate_ell_coeffs {
    bn183_Fq2 ell_0;
    bn183_Fq2 ell_VW;
    bn183_Fq2 ell_VV;

    bool operator==(const bn183_ate_ell_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn183_ate_ell_coeffs &c);
    friend std::istream& operator>>(std::istream &in, bn183_ate_ell_coeffs &c);
};

struct bn183_ate_G2_precomp {
    bn183_Fq2 QX;
    bn183_Fq2 QY;
    std::vector<bn183_ate_ell_coeffs> coeffs;

    bool operator==(const bn183_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn183_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, bn183_ate_G2_precomp &prec_Q);
};

bn183_ate_G1_precomp bn183_ate_precompute_G1(const bn183_G1& P);
bn183_ate_G2_precomp bn183_ate_precompute_G2(const bn183_G2& Q);

bn183_Fq12 bn183_ate_miller_loop(const bn183_ate_G1_precomp &prec_P,
                              const bn183_ate_G2_precomp &prec_Q);
bn183_Fq12 bn183_ate_double_miller_loop(const bn183_ate_G1_precomp &prec_P1,
                                     const bn183_ate_G2_precomp &prec_Q1,
                                     const bn183_ate_G1_precomp &prec_P2,
                                     const bn183_ate_G2_precomp &prec_Q2);

bn183_Fq12 bn183_ate_pairing(const bn183_G1& P,
                          const bn183_G2 &Q);
bn183_GT bn183_ate_reduced_pairing(const bn183_G1 &P,
                                 const bn183_G2 &Q);

/* choice of pairing */

typedef bn183_ate_G1_precomp bn183_G1_precomp;
typedef bn183_ate_G2_precomp bn183_G2_precomp;

bn183_G1_precomp bn183_precompute_G1(const bn183_G1& P);

bn183_G2_precomp bn183_precompute_G2(const bn183_G2& Q);

bn183_Fq12 bn183_miller_loop(const bn183_G1_precomp &prec_P,
                          const bn183_G2_precomp &prec_Q);

bn183_Fq12 bn183_double_miller_loop(const bn183_G1_precomp &prec_P1,
                                 const bn183_G2_precomp &prec_Q1,
                                 const bn183_G1_precomp &prec_P2,
                                 const bn183_G2_precomp &prec_Q2);

bn183_Fq12 bn183_pairing(const bn183_G1& P,
                      const bn183_G2 &Q);

bn183_GT bn183_reduced_pairing(const bn183_G1 &P,
                             const bn183_G2 &Q);

bn183_GT bn183_affine_reduced_pairing(const bn183_G1 &P,
                                    const bn183_G2 &Q);

} // namespace libff
#endif // BN183_PAIRING_HPP_
