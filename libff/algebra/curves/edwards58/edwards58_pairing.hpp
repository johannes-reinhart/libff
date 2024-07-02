#ifndef EDWARDS58_PAIRING_HPP_
#define EDWARDS58_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/edwards58/edwards58_init.hpp>

namespace libff {

/* final exponentiation */

edwards58_Fq6 edwards58_final_exponentiation_last_chunk(const edwards58_Fq6 &elt,
                                                    const edwards58_Fq6 &elt_inv);
edwards58_Fq6 edwards58_final_exponentiation_first_chunk(const edwards58_Fq6 &elt,
                                                     const edwards58_Fq6 &elt_inv);
edwards58_GT edwards58_final_exponentiation(const edwards58_Fq6 &elt);

/* Tate pairing */

struct edwards58_Fq_conic_coefficients {
    edwards58_Fq c_ZZ;
    edwards58_Fq c_XY;
    edwards58_Fq c_XZ;

    bool operator==(const edwards58_Fq_conic_coefficients &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards58_Fq_conic_coefficients &cc);
    friend std::istream& operator>>(std::istream &in, edwards58_Fq_conic_coefficients &cc);
};
typedef std::vector<edwards58_Fq_conic_coefficients> edwards58_tate_G1_precomp;

std::ostream& operator<<(std::ostream& out, const edwards58_tate_G1_precomp &prec_P);
std::istream& operator>>(std::istream& in, edwards58_tate_G1_precomp &prec_P);

struct edwards58_tate_G2_precomp {
    edwards58_Fq3 y0, eta;

    bool operator==(const edwards58_tate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards58_tate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, edwards58_tate_G2_precomp &prec_Q);
};

edwards58_tate_G1_precomp edwards58_tate_precompute_G1(const edwards58_G1& P);
edwards58_tate_G2_precomp edwards58_tate_precompute_G2(const edwards58_G2& Q);

edwards58_Fq6 edwards58_tate_miller_loop(const edwards58_tate_G1_precomp &prec_P,
                                     const edwards58_tate_G2_precomp &prec_Q);

edwards58_Fq6 edwards58_tate_pairing(const edwards58_G1& P,
                                 const edwards58_G2 &Q);
edwards58_GT edwards58_tate_reduced_pairing(const edwards58_G1 &P,
                                        const edwards58_G2 &Q);

/* ate pairing */

struct edwards58_Fq3_conic_coefficients {
    edwards58_Fq3 c_ZZ;
    edwards58_Fq3 c_XY;
    edwards58_Fq3 c_XZ;

    bool operator==(const edwards58_Fq3_conic_coefficients &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards58_Fq3_conic_coefficients &cc);
    friend std::istream& operator>>(std::istream &in, edwards58_Fq3_conic_coefficients &cc);
};
typedef std::vector<edwards58_Fq3_conic_coefficients> edwards58_ate_G2_precomp;

std::ostream& operator<<(std::ostream& out, const edwards58_ate_G2_precomp &prec_Q);
std::istream& operator>>(std::istream& in, edwards58_ate_G2_precomp &prec_Q);

struct edwards58_ate_G1_precomp {
    edwards58_Fq P_XY;
    edwards58_Fq P_XZ;
    edwards58_Fq P_ZZplusYZ;

    bool operator==(const edwards58_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards58_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, edwards58_ate_G1_precomp &prec_P);
};

edwards58_ate_G1_precomp edwards58_ate_precompute_G1(const edwards58_G1& P);
edwards58_ate_G2_precomp edwards58_ate_precompute_G2(const edwards58_G2& Q);

edwards58_Fq6 edwards58_ate_miller_loop(const edwards58_ate_G1_precomp &prec_P,
                                    const edwards58_ate_G2_precomp &prec_Q);
edwards58_Fq6 edwards58_ate_double_miller_loop(const edwards58_ate_G1_precomp &prec_P1,
                                           const edwards58_ate_G2_precomp &prec_Q1,
                                           const edwards58_ate_G1_precomp &prec_P2,
                                           const edwards58_ate_G2_precomp &prec_Q2);

edwards58_Fq6 edwards58_ate_pairing(const edwards58_G1& P,
                                const edwards58_G2 &Q);
edwards58_GT edwards58_ate_reduced_pairing(const edwards58_G1 &P,
                                       const edwards58_G2 &Q);

/* choice of pairing */

typedef edwards58_ate_G1_precomp edwards58_G1_precomp;
typedef edwards58_ate_G2_precomp edwards58_G2_precomp;

edwards58_G1_precomp edwards58_precompute_G1(const edwards58_G1& P);
edwards58_G2_precomp edwards58_precompute_G2(const edwards58_G2& Q);

edwards58_Fq6 edwards58_miller_loop(const edwards58_G1_precomp &prec_P,
                                const edwards58_G2_precomp &prec_Q);

edwards58_Fq6 edwards58_double_miller_loop(const edwards58_G1_precomp &prec_P1,
                                       const edwards58_G2_precomp &prec_Q1,
                                       const edwards58_G1_precomp &prec_P2,
                                       const edwards58_G2_precomp &prec_Q2);

edwards58_Fq6 edwards58_pairing(const edwards58_G1& P,
                            const edwards58_G2 &Q);

edwards58_GT edwards58_reduced_pairing(const edwards58_G1 &P,
                                   const edwards58_G2 &Q);

} // namespace libff
#endif // EDWARDS58_PAIRING_HPP_
