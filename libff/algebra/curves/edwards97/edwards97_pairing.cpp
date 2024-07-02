/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <cassert>

#include <libff/algebra/curves/edwards97/edwards97_g1.hpp>
#include <libff/algebra/curves/edwards97/edwards97_g2.hpp>
#include <libff/algebra/curves/edwards97/edwards97_init.hpp>
#include <libff/algebra/curves/edwards97/edwards97_pairing.hpp>
#include <libff/common/profiling.hpp>

namespace libff {

using std::size_t;

bool edwards97_Fq_conic_coefficients::operator==(const edwards97_Fq_conic_coefficients &other) const
{
    return (this->c_ZZ == other.c_ZZ &&
            this->c_XY == other.c_XY &&
            this->c_XZ == other.c_XZ);
}

std::ostream& operator<<(std::ostream &out, const edwards97_Fq_conic_coefficients &cc)
{
    out << cc.c_ZZ << OUTPUT_SEPARATOR << cc.c_XY << OUTPUT_SEPARATOR << cc.c_XZ;
    return out;
}

std::istream& operator>>(std::istream &in, edwards97_Fq_conic_coefficients &cc)
{
    in >> cc.c_ZZ;
    consume_OUTPUT_SEPARATOR(in);
    in >> cc.c_XY;
    consume_OUTPUT_SEPARATOR(in);
    in >> cc.c_XZ;
    return in;
}

std::ostream& operator<<(std::ostream& out, const edwards97_tate_G1_precomp &prec_P)
{
    out << prec_P.size() << "\n";
    for (const edwards97_Fq_conic_coefficients &cc : prec_P)
    {
        out << cc << OUTPUT_NEWLINE;
    }

    return out;
}

std::istream& operator>>(std::istream& in, edwards97_tate_G1_precomp &prec_P)
{
    prec_P.clear();

    size_t s;
    in >> s;

    consume_newline(in);
    prec_P.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        edwards97_Fq_conic_coefficients cc;
        in >> cc;
        consume_OUTPUT_NEWLINE(in);
        prec_P.emplace_back(cc);
    }

    return in;
}

bool edwards97_tate_G2_precomp::operator==(const edwards97_tate_G2_precomp &other) const
{
    return (this->y0 == other.y0 &&
            this->eta == other.eta);
}

std::ostream& operator<<(std::ostream &out, const edwards97_tate_G2_precomp &prec_Q)
{
    out << prec_Q.y0 << OUTPUT_SEPARATOR << prec_Q.eta;
    return out;
}

std::istream& operator>>(std::istream &in, edwards97_tate_G2_precomp &prec_Q)
{
    in >> prec_Q.y0;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.eta;
    return in;
}

bool edwards97_Fq3_conic_coefficients::operator==(const edwards97_Fq3_conic_coefficients &other) const
{
    return (this->c_ZZ == other.c_ZZ &&
            this->c_XY == other.c_XY &&
            this->c_XZ == other.c_XZ);
}

std::ostream& operator<<(std::ostream &out, const edwards97_Fq3_conic_coefficients &cc)
{
    out << cc.c_ZZ << OUTPUT_SEPARATOR << cc.c_XY << OUTPUT_SEPARATOR << cc.c_XZ;
    return out;
}

std::istream& operator>>(std::istream &in, edwards97_Fq3_conic_coefficients &cc)
{
    in >> cc.c_ZZ;
    consume_OUTPUT_SEPARATOR(in);
    in >> cc.c_XY;
    consume_OUTPUT_SEPARATOR(in);
    in >> cc.c_XZ;
    return in;
}

std::ostream& operator<<(std::ostream& out, const edwards97_ate_G2_precomp &prec_Q)
{
    out << prec_Q.size() << "\n";
    for (const edwards97_Fq3_conic_coefficients &cc : prec_Q)
    {
        out << cc << OUTPUT_NEWLINE;
    }

    return out;
}

std::istream& operator>>(std::istream& in, edwards97_ate_G2_precomp &prec_Q)
{
    prec_Q.clear();

    size_t s;
    in >> s;

    consume_newline(in);

    prec_Q.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        edwards97_Fq3_conic_coefficients cc;
        in >> cc;
        consume_OUTPUT_NEWLINE(in);
        prec_Q.emplace_back(cc);
    }

    return in;
}

bool edwards97_ate_G1_precomp::operator==(const edwards97_ate_G1_precomp &other) const
{
    return (this->P_XY == other.P_XY &&
            this->P_XZ == other.P_XZ &&
            this->P_ZZplusYZ == other.P_ZZplusYZ);
}

std::ostream& operator<<(std::ostream &out, const edwards97_ate_G1_precomp &prec_P)
{
    out << prec_P.P_XY << OUTPUT_SEPARATOR << prec_P.P_XZ << OUTPUT_SEPARATOR << prec_P.P_ZZplusYZ;

    return out;
}

std::istream& operator>>(std::istream &in, edwards97_ate_G1_precomp &prec_P)
{
    in >> prec_P.P_XY >> prec_P.P_XZ >> prec_P.P_ZZplusYZ;

    return in;
}

/* final exponentiations */
edwards97_Fq6 edwards97_final_exponentiation_last_chunk(const edwards97_Fq6 &elt, const edwards97_Fq6 &elt_inv)
{
    enter_block("Call to edwards97_final_exponentiation_last_chunk");
    const edwards97_Fq6 elt_q = elt.Frobenius_map(1);
    edwards97_Fq6 w1_part = elt_q.cyclotomic_exp(edwards97_final_exponent_last_chunk_w1);
    edwards97_Fq6 w0_part;
    if (edwards97_final_exponent_last_chunk_is_w0_neg)
    {
    	w0_part = elt_inv.cyclotomic_exp(edwards97_final_exponent_last_chunk_abs_of_w0);
    } else {
    	w0_part = elt.cyclotomic_exp(edwards97_final_exponent_last_chunk_abs_of_w0);
    }
    edwards97_Fq6 result = w1_part * w0_part;
    leave_block("Call to edwards97_final_exponentiation_last_chunk");

    return result;
}

edwards97_Fq6 edwards97_final_exponentiation_first_chunk(const edwards97_Fq6 &elt, const edwards97_Fq6 &elt_inv)
{
    enter_block("Call to edwards97_final_exponentiation_first_chunk");

    /* (q^3-1)*(q+1) */

    /* elt_q3 = elt^(q^3) */
    const edwards97_Fq6 elt_q3 = elt.Frobenius_map(3);
    /* elt_q3_over_elt = elt^(q^3-1) */
    const edwards97_Fq6 elt_q3_over_elt = elt_q3 * elt_inv;
    /* alpha = elt^((q^3-1) * q) */
    const edwards97_Fq6 alpha = elt_q3_over_elt.Frobenius_map(1);
    /* beta = elt^((q^3-1)*(q+1) */
    const edwards97_Fq6 beta = alpha * elt_q3_over_elt;
    leave_block("Call to edwards97_final_exponentiation_first_chunk");
    return beta;
}

edwards97_GT edwards97_final_exponentiation(const edwards97_Fq6 &elt)
{
    enter_block("Call to edwards97_final_exponentiation");
    const edwards97_Fq6 elt_inv = elt.inverse();
    const edwards97_Fq6 elt_to_first_chunk = edwards97_final_exponentiation_first_chunk(elt, elt_inv);
    const edwards97_Fq6 elt_inv_to_first_chunk = edwards97_final_exponentiation_first_chunk(elt_inv, elt);
    edwards97_GT result = edwards97_final_exponentiation_last_chunk(elt_to_first_chunk, elt_inv_to_first_chunk);
    leave_block("Call to edwards97_final_exponentiation");

    return result;
}

edwards97_tate_G2_precomp edwards97_tate_precompute_G2(const edwards97_G2& Q)
{
    enter_block("Call to edwards97_tate_precompute_G2");
    edwards97_G2 Qcopy = Q;
    Qcopy.to_affine_coordinates();
    edwards97_tate_G2_precomp result;
    result.y0 = Qcopy.Y * Qcopy.Z.inverse(); // Y/Z
    result.eta = (Qcopy.Z+Qcopy.Y) * edwards97_Fq6::mul_by_non_residue(Qcopy.X).inverse(); // (Z+Y)/(nqr*X)
    leave_block("Call to edwards97_tate_precompute_G2");

    return result;
}

struct extended_edwards97_G1_projective {
    edwards97_Fq X;
    edwards97_Fq Y;
    edwards97_Fq Z;
    edwards97_Fq T;

    void print() const
    {
        printf("extended edwards97_G1 projective X/Y/Z/T:\n");
        X.print();
        Y.print();
        Z.print();
        T.print();
    }

    void test_invariant()
    {
        assert(T*Z == X*Y);
    }
};

void doubling_step_for_miller_loop(extended_edwards97_G1_projective &current,
                                   edwards97_Fq_conic_coefficients &cc)
    // From Arene: Faster computation of the Tate pairing
{
    const edwards97_Fq &X = current.X, &Y = current.Y, &Z = current.Z, &T = current.T;
    const edwards97_Fq A = X.squared();     // A    = X1^2
    const edwards97_Fq B = Y.squared();     // B    = Y1^2
    const edwards97_Fq C = Z.squared();     // C    = Z1^2
    const edwards97_Fq D = (X+Y).squared(); // D    = (X1+Y1)^2
    const edwards97_Fq E = (Y+Z).squared(); // E    = (Y1+Z1)^2
    const edwards97_Fq F = D-(A+B);         // F    = D-(A+B)
    const edwards97_Fq G = E-(B+C);         // G    = E-(B+C)
    const edwards97_Fq H = edwards97_G1::coeff_a * A; // H    = a*A
    const edwards97_Fq I = H+B;             // I    = H+B
    const edwards97_Fq J = C-I;             // J    = C-I
    const edwards97_Fq K = J+C;             // K    = J+C

    cc.c_ZZ = Y*(T-X);            // c_ZZ = 2*Y1*(T1-X1)
    cc.c_ZZ = cc.c_ZZ + cc.c_ZZ;

    cc.c_XY = J+J+G;              // c_XY = 2*J+G
    cc.c_XZ = edwards97_G1::coeff_a*X*T-B;              // c_XZ = 2*(aX1*T1-B)
    cc.c_XZ = cc.c_XZ + cc.c_XZ;

    current.X = F*K;              // X3 = F*K
    current.Y = I*(B-H);          // Y3 = I*(B-H)
    current.Z = I*K;              // Z3 = I*K
    current.T = F*(B-H);          // T3 = F*(B-H)

#ifdef DEBUG
    current.test_invariant();
#endif
}

void full_addition_step_for_miller_loop(const extended_edwards97_G1_projective &base,
                                        extended_edwards97_G1_projective &current,
                                        edwards97_Fq_conic_coefficients &cc)
{
    const edwards97_Fq &X1 = current.X, &Y1 = current.Y, &Z1 = current.Z, &T1 = current.T;
    const edwards97_Fq &X2 = base.X, &Y2 =  base.Y, &Z2 = base.Z, &T2 = base.T;

    const edwards97_Fq A = X1*X2;               // A    = X1*X2
    const edwards97_Fq B = Y1*Y2;               // B    = Y1*Y2
    const edwards97_Fq C = Z1*T2;               // C    = Z1*T2
    const edwards97_Fq D = T1*Z2;               // D    = T1*Z2
    const edwards97_Fq E = D+C;                 // E    = D+C
    const edwards97_Fq F = (X1-Y1)*(X2+Y2)+B-A; // F    = (X1-Y1)*(X2+Y2)+B-A
    const edwards97_Fq G = B + edwards97_G1::coeff_a*A;               // G    = B + aA
    const edwards97_Fq H = D-C;                 // H    = D-C
    const edwards97_Fq I = T1*T2;               // I    = T1*T2

    cc.c_ZZ = (T1-X1)*(T2+X2)-I+A;    // c_ZZ = (T1-X1)*(T2+X2)-I+A
    cc.c_XY = X1*Z2-X2*Z1+F;          // c_XY = X1*Z2-X2*Z1+F
    cc.c_XZ = (Y1-T1)*(Y2+T2)-B+I-H;  // c_XZ = (Y1-T1)*(Y2+T2)-B+I-H
    current.X = E*F;                  // X3   = E*F
    current.Y = G*H;                  // Y3   = G*H
    current.Z = F*G;                  // Z3   = F*G
    current.T = E*H;                  // T3   = E*H

#ifdef DEBUG
    current.test_invariant();
#endif
}

void mixed_addition_step_for_miller_loop(const extended_edwards97_G1_projective &base,
                                         extended_edwards97_G1_projective &current,
                                         edwards97_Fq_conic_coefficients &cc)
    //http://www.hyperelliptic.org/EFD/g1p/auto-twisted-extended.html#addition-madd-2008-hwcd-2
{
    const edwards97_Fq &X1 = current.X, &Y1 = current.Y, &Z1 = current.Z, &T1 = current.T;
    const edwards97_Fq &X2 = base.X, &Y2 =  base.Y, &T2 = base.T;

    const edwards97_Fq A = X1*X2;               // A    = X1*X2
    const edwards97_Fq B = Y1*Y2;               // B    = Y1*Y2
    const edwards97_Fq C = Z1*T2;               // C    = Z1*T2
    const edwards97_Fq D = T1;                  // D    = T1*Z2
    const edwards97_Fq E = D+C;                 // E    = D+C
    const edwards97_Fq F = (X1-Y1)*(X2+Y2)+B-A; // F    = (X1-Y1)*(X2+Y2)+B-A
    const edwards97_Fq G = B + edwards97_G1::coeff_a*A;               // G    = B + a*A
    const edwards97_Fq H = D-C;                 // H    = D-C
    const edwards97_Fq I = T1*T2;               // I    = T1*T2

    cc.c_ZZ = (T1-X1)*(T2+X2)-I+A;    // c_ZZ = (T1-X1)*(T2+X2)-I+A
    cc.c_XY = X1-X2*Z1+F;             // c_XY = X1*Z2-X2*Z1+F
    cc.c_XZ = (Y1-T1)*(Y2+T2)-B+I-H;  // c_XZ = (Y1-T1)*(Y2+T2)-B+I-H
    current.X = E*F;                  // X3   = E*F
    current.Y = G*H;                  // Y3   = G*H
    current.Z = F*G;                  // Z3   = F*G
    current.T = E*H;                  // T3   = E*H

#ifdef DEBUG
    current.test_invariant();
#endif
}

edwards97_tate_G1_precomp edwards97_tate_precompute_G1(const edwards97_G1& P)
{
    enter_block("Call to edwards97_tate_precompute_G1");
    edwards97_tate_G1_precomp result;

    edwards97_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();

    extended_edwards97_G1_projective P_ext;
    P_ext.X = Pcopy.X;
    P_ext.Y = Pcopy.Y;
    P_ext.Z = Pcopy.Z;
    P_ext.T = Pcopy.X*Pcopy.Y;

    extended_edwards97_G1_projective R = P_ext;

    bool found_one = false;
    for (long i = edwards97_modulus_r.max_bits(); i >= 0; --i)
    {
        const bool bit = edwards97_modulus_r.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           edwards97_modulus_r (skipping leading zeros) in MSB to LSB
           order */
        edwards97_Fq_conic_coefficients cc;
        doubling_step_for_miller_loop(R, cc);
        result.push_back(cc);

        if (bit)
        {
            mixed_addition_step_for_miller_loop(P_ext, R, cc);
            result.push_back(cc);
        }
    }

    leave_block("Call to edwards97_tate_precompute_G1");
    return result;
}

edwards97_Fq6 edwards97_tate_miller_loop(const edwards97_tate_G1_precomp &prec_P,
                          const edwards97_tate_G2_precomp &prec_Q)
{
    enter_block("Call to edwards97_tate_miller_loop");

    edwards97_Fq6 f = edwards97_Fq6::one();

    bool found_one = false;
    size_t idx = 0;
    for (long i = (long) edwards97_modulus_r.max_bits()-1; i >= 0; --i)
    {
        const bool bit = edwards97_modulus_r.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           edwards97_modulus_r (skipping leading zeros) in MSB to LSB
           order */
        edwards97_Fq_conic_coefficients cc = prec_P[idx++];
        edwards97_Fq6 g_RR_at_Q = edwards97_Fq6(edwards97_Fq3(cc.c_XZ, edwards97_Fq(0L), edwards97_Fq(0L)) + cc.c_XY * prec_Q.y0,
                                            cc.c_ZZ * prec_Q.eta);
        f = f.squared() * g_RR_at_Q;
        if (bit)
        {
            cc = prec_P[idx++];

            edwards97_Fq6 g_RP_at_Q = edwards97_Fq6(edwards97_Fq3(cc.c_XZ, edwards97_Fq(0L), edwards97_Fq(0L)) + cc.c_XY * prec_Q.y0,
                                                cc.c_ZZ * prec_Q.eta);
            f = f * g_RP_at_Q;
        }
    }
    leave_block("Call to edwards97_tate_miller_loop");

    return f;
}

edwards97_Fq6 edwards97_tate_pairing(const edwards97_G1& P, const edwards97_G2 &Q)
{
    enter_block("Call to edwards97_tate_pairing");
    edwards97_tate_G1_precomp prec_P = edwards97_tate_precompute_G1(P);
    edwards97_tate_G2_precomp prec_Q = edwards97_tate_precompute_G2(Q);
    edwards97_Fq6 result = edwards97_tate_miller_loop(prec_P, prec_Q);
    leave_block("Call to edwards97_tate_pairing");
    return result;
}

edwards97_GT edwards97_tate_reduced_pairing(const edwards97_G1 &P, const edwards97_G2 &Q)
{
    enter_block("Call to edwards97_tate_reduced_pairing");
    const edwards97_Fq6 f = edwards97_tate_pairing(P, Q);
    const edwards97_GT result = edwards97_final_exponentiation(f);
    leave_block("Call to edwards97_tate_reduce_pairing");
    return result;
}

struct extended_edwards97_G2_projective {
    edwards97_Fq3 X;
    edwards97_Fq3 Y;
    edwards97_Fq3 Z;
    edwards97_Fq3 T;

    void print() const
        {
            printf("extended edwards97_G2 projective X/Y/Z/T:\n");
            X.print();
            Y.print();
            Z.print();
            T.print();
        }

    void test_invariant()
        {
            assert(T*Z == X*Y);
        }
};

void doubling_step_for_flipped_miller_loop(extended_edwards97_G2_projective &current,
                                           edwards97_Fq3_conic_coefficients &cc)
{
    const edwards97_Fq3 &X = current.X, &Y = current.Y, &Z = current.Z, &T = current.T;
    const edwards97_Fq3 A = X.squared();     // A    = X1^2
    const edwards97_Fq3 B = Y.squared();     // B    = Y1^2
    const edwards97_Fq3 C = Z.squared();     // C    = Z1^2
    const edwards97_Fq3 D = (X+Y).squared(); // D    = (X1+Y1)^2
    const edwards97_Fq3 E = (Y+Z).squared(); // E    = (Y1+Z1)^2
    const edwards97_Fq3 F = D-(A+B);         // F    = D-(A+B)
    const edwards97_Fq3 G = E-(B+C);         // G    = E-(B+C)
    const edwards97_Fq3 H = edwards97_G2::mul_by_a(A);  // H    = twisted_a * A
    const edwards97_Fq3 I = H+B;             // I    = H+B
    const edwards97_Fq3 J = C-I;             // J    = C-I
    const edwards97_Fq3 K = J+C;             // K    = J+C

    cc.c_ZZ = Y*(T-X);            // c_ZZ = 2*Y1*(T1-X1)
    cc.c_ZZ = cc.c_ZZ + cc.c_ZZ;

    // c_XY = 2*(C-edwards_a * A * delta_3-B)+G
    // cc.c_XY = C - edwards97_G2::mul_by_a(A) - B;
    cc.c_XY = J + J + G;

    // c_XZ = 2*(edwards_a*X1*T1*delta_3-B)
    cc.c_XZ = edwards97_G2::mul_by_a(X * T) - B;
    cc.c_XZ = cc.c_XZ + cc.c_XZ;

    current.X = F*K;              // X3 = F*K
    current.Y = I*(B-H);          // Y3 = I*(B-H)
    current.Z = I*K;              // Z3 = I*K
    current.T = F*(B-H);          // T3 = F*(B-H)
#ifdef DEBUG
    current.test_invariant();
#endif
}

void full_addition_step_for_flipped_miller_loop(const extended_edwards97_G2_projective &base,
                                                extended_edwards97_G2_projective &current,
                                                edwards97_Fq3_conic_coefficients &cc)
{
    const edwards97_Fq3 &X1 = current.X, &Y1 = current.Y, &Z1 = current.Z, &T1 = current.T;
    const edwards97_Fq3 &X2 = base.X, &Y2 =  base.Y, &Z2 = base.Z, &T2 = base.T;

    const edwards97_Fq3 A = X1*X2;               // A    = X1*X2
    const edwards97_Fq3 B = Y1*Y2;               // B    = Y1*Y2
    const edwards97_Fq3 C = Z1*T2;               // C    = Z1*T2
    const edwards97_Fq3 D = T1*Z2;               // D    = T1*Z2
    const edwards97_Fq3 E = D+C;                 // E    = D+C
    const edwards97_Fq3 F = (X1-Y1)*(X2+Y2)+B-A; // F    = (X1-Y1)*(X2+Y2)+B-A
    // G = B + twisted_edwards97_a * A
    const edwards97_Fq3 G = B + edwards97_G2::mul_by_a(A);
    const edwards97_Fq3 H = D-C;                 // H    = D-C
    const edwards97_Fq3 I = T1*T2;               // I    = T1*T2

    // c_ZZ = delta_3* ((T1-X1)*(T2+X2)-I+A)
    cc.c_ZZ = edwards97_G2::mul_by_twist((T1-X1)*(T2+X2)-I+A);

    cc.c_XY = X1*Z2-X2*Z1+F;          // c_XY = X1*Z2-X2*Z1+F
    cc.c_XZ = (Y1-T1)*(Y2+T2)-B+I-H;  // c_XZ = (Y1-T1)*(Y2+T2)-B+I-H
    current.X = E*F;                  // X3   = E*F
    current.Y = G*H;                  // Y3   = G*H
    current.Z = F*G;                  // Z3   = F*G
    current.T = E*H;                  // T3   = E*H

#ifdef DEBUG
    current.test_invariant();
#endif
}

void mixed_addition_step_for_flipped_miller_loop(const extended_edwards97_G2_projective &base,
                                                 extended_edwards97_G2_projective &current,
                                                 edwards97_Fq3_conic_coefficients &cc)
{
    const edwards97_Fq3 &X1 = current.X, &Y1 = current.Y, &Z1 = current.Z, &T1 = current.T;
    const edwards97_Fq3 &X2 = base.X, &Y2 =  base.Y, &T2 = base.T;

    const edwards97_Fq3 A = X1*X2;               // A    = X1*X2
    const edwards97_Fq3 B = Y1*Y2;               // B    = Y1*Y2
    const edwards97_Fq3 C = Z1*T2;               // C    = Z1*T2
    const edwards97_Fq3 E = T1+C;                // E    = T1+C
    const edwards97_Fq3 F = (X1-Y1)*(X2+Y2)+B-A; // F    = (X1-Y1)*(X2+Y2)+B-A
    // G = B + twisted_edwards97_a * A
    const edwards97_Fq3 G = B + edwards97_G2::mul_by_a(A);
    const edwards97_Fq3 H = T1-C;                // H    = T1-C
    const edwards97_Fq3 I = T1*T2;               // I    = T1*T2

    // c_ZZ = delta_3* ((T1-X1)*(T2+X2)-I+A)
    cc.c_ZZ = edwards97_G2::mul_by_twist((T1-X1)*(T2+X2)-I+A);

    cc.c_XY = X1-X2*Z1+F;             // c_XY = X1*Z2-X2*Z1+F
    cc.c_XZ = (Y1-T1)*(Y2+T2)-B+I-H;  // c_XZ = (Y1-T1)*(Y2+T2)-B+I-H
    current.X = E*F;                  // X3   = E*F
    current.Y = G*H;                  // Y3   = G*H
    current.Z = F*G;                  // Z3   = F*G
    current.T = E*H;                  // T3   = E*H

#ifdef DEBUG
    current.test_invariant();
#endif
}

edwards97_ate_G1_precomp edwards97_ate_precompute_G1(const edwards97_G1& P)
{
    enter_block("Call to edwards97_ate_precompute_G1");
    edwards97_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();
    edwards97_ate_G1_precomp result;
    result.P_XY = Pcopy.X*Pcopy.Y;
    result.P_XZ = Pcopy.X; // P.X * P.Z but P.Z = 1
    result.P_ZZplusYZ = (edwards97_Fq::one() + Pcopy.Y); // (P.Z + P.Y) * P.Z but P.Z = 1
    leave_block("Call to edwards97_ate_precompute_G1");
    return result;
}

edwards97_ate_G2_precomp edwards97_ate_precompute_G2(const edwards97_G2& Q)
{
    enter_block("Call to edwards97_ate_precompute_G2");
    const bigint<edwards97_Fr::num_limbs> &loop_count = edwards97_ate_loop_count;
    edwards97_ate_G2_precomp result;

    edwards97_G2 Qcopy(Q);
    Qcopy.to_affine_coordinates();

    extended_edwards97_G2_projective Q_ext;
    Q_ext.X = Qcopy.X;
    Q_ext.Y = Qcopy.Y;
    Q_ext.Z = Qcopy.Z;
    Q_ext.T = Qcopy.X*Qcopy.Y;

    extended_edwards97_G2_projective R = Q_ext;

    bool found_one = false;
    for (long i = (long) loop_count.max_bits()-1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        edwards97_Fq3_conic_coefficients cc;
        doubling_step_for_flipped_miller_loop(R, cc);
        result.push_back(cc);
        if (bit)
        {
            mixed_addition_step_for_flipped_miller_loop(Q_ext, R, cc);
            result.push_back(cc);
        }
    }

    leave_block("Call to edwards97_ate_precompute_G2");
    return result;
}

edwards97_Fq6 edwards97_ate_miller_loop(const edwards97_ate_G1_precomp &prec_P,
                                    const edwards97_ate_G2_precomp &prec_Q)
{
    enter_block("Call to edwards97_ate_miller_loop");
    const bigint<edwards97_Fr::num_limbs> &loop_count = edwards97_ate_loop_count;

    edwards97_Fq6 f = edwards97_Fq6::one();

    bool found_one = false;
    size_t idx = 0;
    for (long i = (long) loop_count.max_bits()-1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           edwards_param_p (skipping leading zeros) in MSB to LSB
           order */
        edwards97_Fq3_conic_coefficients cc = prec_Q[idx++];

        edwards97_Fq6 g_RR_at_P = edwards97_Fq6(prec_P.P_XY * cc.c_XY + prec_P.P_XZ * cc.c_XZ,
                                            prec_P.P_ZZplusYZ * cc.c_ZZ);
        f = f.squared() * g_RR_at_P;
        if (bit)
        {
            cc = prec_Q[idx++];
            edwards97_Fq6 g_RQ_at_P = edwards97_Fq6(prec_P.P_ZZplusYZ * cc.c_ZZ,
                                                prec_P.P_XY * cc.c_XY + prec_P.P_XZ * cc.c_XZ);
            f = f * g_RQ_at_P;
        }
    }
    leave_block("Call to edwards97_ate_miller_loop");

    return f;
}

edwards97_Fq6 edwards97_ate_double_miller_loop(const edwards97_ate_G1_precomp &prec_P1,
                                           const edwards97_ate_G2_precomp &prec_Q1,
                                           const edwards97_ate_G1_precomp &prec_P2,
                                           const edwards97_ate_G2_precomp &prec_Q2)
{
    enter_block("Call to edwards97_ate_double_miller_loop");
    const bigint<edwards97_Fr::num_limbs> &loop_count = edwards97_ate_loop_count;

    edwards97_Fq6 f = edwards97_Fq6::one();

    bool found_one = false;
    size_t idx = 0;
    for (long i = (long) loop_count.max_bits()-1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           edwards_param_p (skipping leading zeros) in MSB to LSB
           order */
        edwards97_Fq3_conic_coefficients cc1 = prec_Q1[idx];
        edwards97_Fq3_conic_coefficients cc2 = prec_Q2[idx];
        ++idx;

        edwards97_Fq6 g_RR_at_P1 = edwards97_Fq6(prec_P1.P_XY * cc1.c_XY + prec_P1.P_XZ * cc1.c_XZ,
                                             prec_P1.P_ZZplusYZ * cc1.c_ZZ);

        edwards97_Fq6 g_RR_at_P2 = edwards97_Fq6(prec_P2.P_XY * cc2.c_XY + prec_P2.P_XZ * cc2.c_XZ,
                                             prec_P2.P_ZZplusYZ * cc2.c_ZZ);
        f = f.squared() * g_RR_at_P1 * g_RR_at_P2;

        if (bit)
        {
            cc1 = prec_Q1[idx];
            cc2 = prec_Q2[idx];
            ++idx;
            edwards97_Fq6 g_RQ_at_P1 = edwards97_Fq6(prec_P1.P_ZZplusYZ * cc1.c_ZZ,
                                                 prec_P1.P_XY * cc1.c_XY + prec_P1.P_XZ * cc1.c_XZ);
            edwards97_Fq6 g_RQ_at_P2 = edwards97_Fq6(prec_P2.P_ZZplusYZ * cc2.c_ZZ,
                                                 prec_P2.P_XY * cc2.c_XY + prec_P2.P_XZ * cc2.c_XZ);
            f = f * g_RQ_at_P1 * g_RQ_at_P2;
        }
    }
    leave_block("Call to edwards97_ate_double_miller_loop");

    return f;
}

edwards97_Fq6 edwards97_ate_pairing(const edwards97_G1& P, const edwards97_G2 &Q)
{
    enter_block("Call to edwards97_ate_pairing");
    edwards97_ate_G1_precomp prec_P = edwards97_ate_precompute_G1(P);
    edwards97_ate_G2_precomp prec_Q = edwards97_ate_precompute_G2(Q);
    edwards97_Fq6 result = edwards97_ate_miller_loop(prec_P, prec_Q);
    leave_block("Call to edwards97_ate_pairing");
    return result;
}

edwards97_GT edwards97_ate_reduced_pairing(const edwards97_G1 &P, const edwards97_G2 &Q)
{
    enter_block("Call to edwards97_ate_reduced_pairing");
    const edwards97_Fq6 f = edwards97_ate_pairing(P, Q);
    const edwards97_GT result = edwards97_final_exponentiation(f);
    leave_block("Call to edwards97_ate_reduced_pairing");
    return result;
}

edwards97_G1_precomp edwards97_precompute_G1(const edwards97_G1& P)
{
    return edwards97_ate_precompute_G1(P);
}

edwards97_G2_precomp edwards97_precompute_G2(const edwards97_G2& Q)
{
    return edwards97_ate_precompute_G2(Q);
}

edwards97_Fq6 edwards97_miller_loop(const edwards97_G1_precomp &prec_P,
                                const edwards97_G2_precomp &prec_Q)
{
    return edwards97_ate_miller_loop(prec_P, prec_Q);
}

edwards97_Fq6 edwards97_double_miller_loop(const edwards97_G1_precomp &prec_P1,
                                       const edwards97_G2_precomp &prec_Q1,
                                       const edwards97_G1_precomp &prec_P2,
                                       const edwards97_G2_precomp &prec_Q2)
{
    return edwards97_ate_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

edwards97_Fq6 edwards97_pairing(const edwards97_G1& P,
                            const edwards97_G2 &Q)
{
    return edwards97_ate_pairing(P, Q);
}

edwards97_GT edwards97_reduced_pairing(const edwards97_G1 &P,
                                   const edwards97_G2 &Q)
{
    return edwards97_ate_reduced_pairing(P, Q);
}
} // namespace libff
