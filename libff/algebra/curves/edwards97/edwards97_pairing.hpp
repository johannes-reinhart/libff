#ifndef EDWARDS97_PAIRING_HPP_
#define EDWARDS97_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/edwards97/edwards97_init.hpp>

namespace libff {

/* final exponentiation */

edwards97_Fq6 edwards97_final_exponentiation_last_chunk(const edwards97_Fq6 &elt,
                                                    const edwards97_Fq6 &elt_inv);
edwards97_Fq6 edwards97_final_exponentiation_first_chunk(const edwards97_Fq6 &elt,
                                                     const edwards97_Fq6 &elt_inv);
edwards97_GT edwards97_final_exponentiation(const edwards97_Fq6 &elt);

/* Tate pairing */

struct edwards97_Fq_conic_coefficients {
    edwards97_Fq c_ZZ;
    edwards97_Fq c_XY;
    edwards97_Fq c_XZ;

    bool operator==(const edwards97_Fq_conic_coefficients &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards97_Fq_conic_coefficients &cc);
    friend std::istream& operator>>(std::istream &in, edwards97_Fq_conic_coefficients &cc);
};
typedef std::vector<edwards97_Fq_conic_coefficients> edwards97_tate_G1_precomp;

std::ostream& operator<<(std::ostream& out, const edwards97_tate_G1_precomp &prec_P);
std::istream& operator>>(std::istream& in, edwards97_tate_G1_precomp &prec_P);

struct edwards97_tate_G2_precomp {
    edwards97_Fq3 y0, eta;

    bool operator==(const edwards97_tate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards97_tate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, edwards97_tate_G2_precomp &prec_Q);
};

edwards97_tate_G1_precomp edwards97_tate_precompute_G1(const edwards97_G1& P);
edwards97_tate_G2_precomp edwards97_tate_precompute_G2(const edwards97_G2& Q);

edwards97_Fq6 edwards97_tate_miller_loop(const edwards97_tate_G1_precomp &prec_P,
                                     const edwards97_tate_G2_precomp &prec_Q);

edwards97_Fq6 edwards97_tate_pairing(const edwards97_G1& P,
                                 const edwards97_G2 &Q);
edwards97_GT edwards97_tate_reduced_pairing(const edwards97_G1 &P,
                                        const edwards97_G2 &Q);

/* ate pairing */

struct edwards97_Fq3_conic_coefficients {
    edwards97_Fq3 c_ZZ;
    edwards97_Fq3 c_XY;
    edwards97_Fq3 c_XZ;

    bool operator==(const edwards97_Fq3_conic_coefficients &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards97_Fq3_conic_coefficients &cc);
    friend std::istream& operator>>(std::istream &in, edwards97_Fq3_conic_coefficients &cc);
};
typedef std::vector<edwards97_Fq3_conic_coefficients> edwards97_ate_G2_precomp;

std::ostream& operator<<(std::ostream& out, const edwards97_ate_G2_precomp &prec_Q);
std::istream& operator>>(std::istream& in, edwards97_ate_G2_precomp &prec_Q);

struct edwards97_ate_G1_precomp {
    edwards97_Fq P_XY;
    edwards97_Fq P_XZ;
    edwards97_Fq P_ZZplusYZ;

    bool operator==(const edwards97_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards97_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, edwards97_ate_G1_precomp &prec_P);
};

edwards97_ate_G1_precomp edwards97_ate_precompute_G1(const edwards97_G1& P);
edwards97_ate_G2_precomp edwards97_ate_precompute_G2(const edwards97_G2& Q);

edwards97_Fq6 edwards97_ate_miller_loop(const edwards97_ate_G1_precomp &prec_P,
                                    const edwards97_ate_G2_precomp &prec_Q);
edwards97_Fq6 edwards97_ate_double_miller_loop(const edwards97_ate_G1_precomp &prec_P1,
                                           const edwards97_ate_G2_precomp &prec_Q1,
                                           const edwards97_ate_G1_precomp &prec_P2,
                                           const edwards97_ate_G2_precomp &prec_Q2);

edwards97_Fq6 edwards97_ate_pairing(const edwards97_G1& P,
                                const edwards97_G2 &Q);
edwards97_GT edwards97_ate_reduced_pairing(const edwards97_G1 &P,
                                       const edwards97_G2 &Q);

/* choice of pairing */

typedef edwards97_ate_G1_precomp edwards97_G1_precomp;
typedef edwards97_ate_G2_precomp edwards97_G2_precomp;

edwards97_G1_precomp edwards97_precompute_G1(const edwards97_G1& P);
edwards97_G2_precomp edwards97_precompute_G2(const edwards97_G2& Q);

edwards97_Fq6 edwards97_miller_loop(const edwards97_G1_precomp &prec_P,
                                const edwards97_G2_precomp &prec_Q);

edwards97_Fq6 edwards97_double_miller_loop(const edwards97_G1_precomp &prec_P1,
                                       const edwards97_G2_precomp &prec_Q1,
                                       const edwards97_G1_precomp &prec_P2,
                                       const edwards97_G2_precomp &prec_Q2);

edwards97_Fq6 edwards97_pairing(const edwards97_G1& P,
                            const edwards97_G2 &Q);

edwards97_GT edwards97_reduced_pairing(const edwards97_G1 &P,
                                   const edwards97_G2 &Q);

} // namespace libff
#endif // EDWARDS97_PAIRING_HPP_
