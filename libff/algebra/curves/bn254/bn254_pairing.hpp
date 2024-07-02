#ifndef BN254_PAIRING_HPP_
#define BN254_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/bn254/bn254_init.hpp>

namespace libff {

/* final exponentiation */

bn254_GT bn254_final_exponentiation(const bn254_Fq12 &elt);

/* ate pairing */

struct bn254_ate_G1_precomp {
    bn254_Fq PX;
    bn254_Fq PY;

    bool operator==(const bn254_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn254_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, bn254_ate_G1_precomp &prec_P);
};

struct bn254_ate_ell_coeffs {
    bn254_Fq2 ell_0;
    bn254_Fq2 ell_VW;
    bn254_Fq2 ell_VV;

    bool operator==(const bn254_ate_ell_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn254_ate_ell_coeffs &c);
    friend std::istream& operator>>(std::istream &in, bn254_ate_ell_coeffs &c);
};

struct bn254_ate_G2_precomp {
    bn254_Fq2 QX;
    bn254_Fq2 QY;
    std::vector<bn254_ate_ell_coeffs> coeffs;

    bool operator==(const bn254_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bn254_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, bn254_ate_G2_precomp &prec_Q);
};

bn254_ate_G1_precomp bn254_ate_precompute_G1(const bn254_G1& P);
bn254_ate_G2_precomp bn254_ate_precompute_G2(const bn254_G2& Q);

bn254_Fq12 bn254_ate_miller_loop(const bn254_ate_G1_precomp &prec_P,
                              const bn254_ate_G2_precomp &prec_Q);
bn254_Fq12 bn254_ate_double_miller_loop(const bn254_ate_G1_precomp &prec_P1,
                                     const bn254_ate_G2_precomp &prec_Q1,
                                     const bn254_ate_G1_precomp &prec_P2,
                                     const bn254_ate_G2_precomp &prec_Q2);

bn254_Fq12 bn254_ate_pairing(const bn254_G1& P,
                          const bn254_G2 &Q);
bn254_GT bn254_ate_reduced_pairing(const bn254_G1 &P,
                                 const bn254_G2 &Q);

/* choice of pairing */

typedef bn254_ate_G1_precomp bn254_G1_precomp;
typedef bn254_ate_G2_precomp bn254_G2_precomp;

bn254_G1_precomp bn254_precompute_G1(const bn254_G1& P);

bn254_G2_precomp bn254_precompute_G2(const bn254_G2& Q);

bn254_Fq12 bn254_miller_loop(const bn254_G1_precomp &prec_P,
                          const bn254_G2_precomp &prec_Q);

bn254_Fq12 bn254_double_miller_loop(const bn254_G1_precomp &prec_P1,
                                 const bn254_G2_precomp &prec_Q1,
                                 const bn254_G1_precomp &prec_P2,
                                 const bn254_G2_precomp &prec_Q2);

bn254_Fq12 bn254_pairing(const bn254_G1& P,
                      const bn254_G2 &Q);

bn254_GT bn254_reduced_pairing(const bn254_G1 &P,
                             const bn254_G2 &Q);

bn254_GT bn254_affine_reduced_pairing(const bn254_G1 &P,
                                    const bn254_G2 &Q);

} // namespace libff
#endif // BN254_PAIRING_HPP_
