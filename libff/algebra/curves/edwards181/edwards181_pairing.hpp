#ifndef EDWARDS181_PAIRING_HPP_
#define EDWARDS181_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/edwards181/edwards181_init.hpp>

namespace libff {

/* final exponentiation */

edwards181_Fq6 edwards181_final_exponentiation_last_chunk(const edwards181_Fq6 &elt,
                                                    const edwards181_Fq6 &elt_inv);
edwards181_Fq6 edwards181_final_exponentiation_first_chunk(const edwards181_Fq6 &elt,
                                                     const edwards181_Fq6 &elt_inv);
edwards181_GT edwards181_final_exponentiation(const edwards181_Fq6 &elt);

/* Tate pairing */

struct edwards181_Fq_conic_coefficients {
    edwards181_Fq c_ZZ;
    edwards181_Fq c_XY;
    edwards181_Fq c_XZ;

    bool operator==(const edwards181_Fq_conic_coefficients &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards181_Fq_conic_coefficients &cc);
    friend std::istream& operator>>(std::istream &in, edwards181_Fq_conic_coefficients &cc);
};
typedef std::vector<edwards181_Fq_conic_coefficients> edwards181_tate_G1_precomp;

std::ostream& operator<<(std::ostream& out, const edwards181_tate_G1_precomp &prec_P);
std::istream& operator>>(std::istream& in, edwards181_tate_G1_precomp &prec_P);

struct edwards181_tate_G2_precomp {
    edwards181_Fq3 y0, eta;

    bool operator==(const edwards181_tate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards181_tate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, edwards181_tate_G2_precomp &prec_Q);
};

edwards181_tate_G1_precomp edwards181_tate_precompute_G1(const edwards181_G1& P);
edwards181_tate_G2_precomp edwards181_tate_precompute_G2(const edwards181_G2& Q);

edwards181_Fq6 edwards181_tate_miller_loop(const edwards181_tate_G1_precomp &prec_P,
                                     const edwards181_tate_G2_precomp &prec_Q);

edwards181_Fq6 edwards181_tate_pairing(const edwards181_G1& P,
                                 const edwards181_G2 &Q);
edwards181_GT edwards181_tate_reduced_pairing(const edwards181_G1 &P,
                                        const edwards181_G2 &Q);

/* ate pairing */

struct edwards181_Fq3_conic_coefficients {
    edwards181_Fq3 c_ZZ;
    edwards181_Fq3 c_XY;
    edwards181_Fq3 c_XZ;

    bool operator==(const edwards181_Fq3_conic_coefficients &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards181_Fq3_conic_coefficients &cc);
    friend std::istream& operator>>(std::istream &in, edwards181_Fq3_conic_coefficients &cc);
};
typedef std::vector<edwards181_Fq3_conic_coefficients> edwards181_ate_G2_precomp;

std::ostream& operator<<(std::ostream& out, const edwards181_ate_G2_precomp &prec_Q);
std::istream& operator>>(std::istream& in, edwards181_ate_G2_precomp &prec_Q);

struct edwards181_ate_G1_precomp {
    edwards181_Fq P_XY;
    edwards181_Fq P_XZ;
    edwards181_Fq P_ZZplusYZ;

    bool operator==(const edwards181_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const edwards181_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, edwards181_ate_G1_precomp &prec_P);
};

edwards181_ate_G1_precomp edwards181_ate_precompute_G1(const edwards181_G1& P);
edwards181_ate_G2_precomp edwards181_ate_precompute_G2(const edwards181_G2& Q);

edwards181_Fq6 edwards181_ate_miller_loop(const edwards181_ate_G1_precomp &prec_P,
                                    const edwards181_ate_G2_precomp &prec_Q);
edwards181_Fq6 edwards181_ate_double_miller_loop(const edwards181_ate_G1_precomp &prec_P1,
                                           const edwards181_ate_G2_precomp &prec_Q1,
                                           const edwards181_ate_G1_precomp &prec_P2,
                                           const edwards181_ate_G2_precomp &prec_Q2);

edwards181_Fq6 edwards181_ate_pairing(const edwards181_G1& P,
                                const edwards181_G2 &Q);
edwards181_GT edwards181_ate_reduced_pairing(const edwards181_G1 &P,
                                       const edwards181_G2 &Q);

/* choice of pairing */

typedef edwards181_ate_G1_precomp edwards181_G1_precomp;
typedef edwards181_ate_G2_precomp edwards181_G2_precomp;

edwards181_G1_precomp edwards181_precompute_G1(const edwards181_G1& P);
edwards181_G2_precomp edwards181_precompute_G2(const edwards181_G2& Q);

edwards181_Fq6 edwards181_miller_loop(const edwards181_G1_precomp &prec_P,
                                const edwards181_G2_precomp &prec_Q);

edwards181_Fq6 edwards181_double_miller_loop(const edwards181_G1_precomp &prec_P1,
                                       const edwards181_G2_precomp &prec_Q1,
                                       const edwards181_G1_precomp &prec_P2,
                                       const edwards181_G2_precomp &prec_Q2);

edwards181_Fq6 edwards181_pairing(const edwards181_G1& P,
                            const edwards181_G2 &Q);

edwards181_GT edwards181_reduced_pairing(const edwards181_G1 &P,
                                   const edwards181_G2 &Q);

} // namespace libff
#endif // EDWARDS181_PAIRING_HPP_
