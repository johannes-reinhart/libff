#ifndef EDWARDS61_PAIRING_HPP_
#define EDWARDS61_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/edwards61/edwards61_init.hpp>

namespace libff {

/* final exponentiation */

edwards61_Fq6 edwards61_final_exponentiation_last_chunk(const edwards61_Fq6 &elt,
                                                    const edwards61_Fq6 &elt_inv);
edwards61_Fq6 edwards61_final_exponentiation_first_chunk(const edwards61_Fq6 &elt,
                                                     const edwards61_Fq6 &elt_inv);
edwards61_GT edwards61_final_exponentiation(const edwards61_Fq6 &elt);

/* Tate pairing */

struct edwards61_Fq_conic_coefficients {
    edwards61_Fq c_ZZ;
    edwards61_Fq c_XY;
    edwards61_Fq c_XZ;

    bool operator==(const edwards61_Fq_conic_coefficients &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards61_Fq_conic_coefficients &cc);
    friend std::istream& operator>>(std::istream &in, edwards61_Fq_conic_coefficients &cc);
};
typedef std::vector<edwards61_Fq_conic_coefficients> edwards61_tate_G1_precomp;

std::ostream& operator<<(std::ostream& out, const edwards61_tate_G1_precomp &prec_P);
std::istream& operator>>(std::istream& in, edwards61_tate_G1_precomp &prec_P);

struct edwards61_tate_G2_precomp {
    edwards61_Fq3 y0, eta;

    bool operator==(const edwards61_tate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards61_tate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, edwards61_tate_G2_precomp &prec_Q);
};

edwards61_tate_G1_precomp edwards61_tate_precompute_G1(const edwards61_G1& P);
edwards61_tate_G2_precomp edwards61_tate_precompute_G2(const edwards61_G2& Q);

edwards61_Fq6 edwards61_tate_miller_loop(const edwards61_tate_G1_precomp &prec_P,
                                     const edwards61_tate_G2_precomp &prec_Q);

edwards61_Fq6 edwards61_tate_pairing(const edwards61_G1& P,
                                 const edwards61_G2 &Q);
edwards61_GT edwards61_tate_reduced_pairing(const edwards61_G1 &P,
                                        const edwards61_G2 &Q);

/* ate pairing */

struct edwards61_Fq3_conic_coefficients {
    edwards61_Fq3 c_ZZ;
    edwards61_Fq3 c_XY;
    edwards61_Fq3 c_XZ;

    bool operator==(const edwards61_Fq3_conic_coefficients &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards61_Fq3_conic_coefficients &cc);
    friend std::istream& operator>>(std::istream &in, edwards61_Fq3_conic_coefficients &cc);
};
typedef std::vector<edwards61_Fq3_conic_coefficients> edwards61_ate_G2_precomp;

std::ostream& operator<<(std::ostream& out, const edwards61_ate_G2_precomp &prec_Q);
std::istream& operator>>(std::istream& in, edwards61_ate_G2_precomp &prec_Q);

struct edwards61_ate_G1_precomp {
    edwards61_Fq P_XY;
    edwards61_Fq P_XZ;
    edwards61_Fq P_ZZplusYZ;

    bool operator==(const edwards61_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards61_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, edwards61_ate_G1_precomp &prec_P);
};

edwards61_ate_G1_precomp edwards61_ate_precompute_G1(const edwards61_G1& P);
edwards61_ate_G2_precomp edwards61_ate_precompute_G2(const edwards61_G2& Q);

edwards61_Fq6 edwards61_ate_miller_loop(const edwards61_ate_G1_precomp &prec_P,
                                    const edwards61_ate_G2_precomp &prec_Q);
edwards61_Fq6 edwards61_ate_double_miller_loop(const edwards61_ate_G1_precomp &prec_P1,
                                           const edwards61_ate_G2_precomp &prec_Q1,
                                           const edwards61_ate_G1_precomp &prec_P2,
                                           const edwards61_ate_G2_precomp &prec_Q2);

edwards61_Fq6 edwards61_ate_pairing(const edwards61_G1& P,
                                const edwards61_G2 &Q);
edwards61_GT edwards61_ate_reduced_pairing(const edwards61_G1 &P,
                                       const edwards61_G2 &Q);

/* choice of pairing */

typedef edwards61_ate_G1_precomp edwards61_G1_precomp;
typedef edwards61_ate_G2_precomp edwards61_G2_precomp;

edwards61_G1_precomp edwards61_precompute_G1(const edwards61_G1& P);
edwards61_G2_precomp edwards61_precompute_G2(const edwards61_G2& Q);

edwards61_Fq6 edwards61_miller_loop(const edwards61_G1_precomp &prec_P,
                                const edwards61_G2_precomp &prec_Q);

edwards61_Fq6 edwards61_double_miller_loop(const edwards61_G1_precomp &prec_P1,
                                       const edwards61_G2_precomp &prec_Q1,
                                       const edwards61_G1_precomp &prec_P2,
                                       const edwards61_G2_precomp &prec_Q2);

edwards61_Fq6 edwards61_pairing(const edwards61_G1& P,
                            const edwards61_G2 &Q);

edwards61_GT edwards61_reduced_pairing(const edwards61_G1 &P,
                                   const edwards61_G2 &Q);

} // namespace libff
#endif // EDWARDS61_PAIRING_HPP_
