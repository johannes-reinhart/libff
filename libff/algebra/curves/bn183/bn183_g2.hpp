#ifndef BN183_G2_HPP_
#define BN183_G2_HPP_
#include <vector>

#include <libff/algebra/curves/bn183/bn183_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class bn183_G2;
std::ostream& operator<<(std::ostream &, const bn183_G2&);
std::istream& operator>>(std::istream &, bn183_G2&);

class bn183_G2 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<std::size_t> wnaf_window_table;
    static std::vector<std::size_t> fixed_base_exp_window_table;
    static bn183_G2 G2_zero;
    static bn183_G2 G2_one;
    static bool initialized;

    typedef bn183_Fq base_field;
    typedef bn183_Fq2 twist_field;
    typedef bn183_Fr scalar_field;

    // Cofactor
    static const mp_size_t h_bitcount = 183;
    static const mp_size_t h_limbs = (h_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
    static bigint<h_limbs> h;

    bn183_Fq2 X, Y, Z;

    // using Jacobian coordinates
    bn183_G2();
    bn183_G2(const bn183_Fq2& X, const bn183_Fq2& Y) : X(X), Y(Y), Z(bn183_Fq2::one()) {};
    bn183_G2(const bn183_Fq2& X, const bn183_Fq2& Y, const bn183_Fq2& Z) : X(X), Y(Y), Z(Z) {};

    static bn183_Fq2 mul_by_b(const bn183_Fq2 &elt);

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const bn183_G2 &other) const;
    bool operator!=(const bn183_G2 &other) const;

    bn183_G2 operator+(const bn183_G2 &other) const;
    bn183_G2 operator-() const;
    bn183_G2 operator-(const bn183_G2 &other) const;

    bn183_G2 add(const bn183_G2 &other) const;
    bn183_G2 mixed_add(const bn183_G2 &other) const;
    bn183_G2 dbl() const;
    bn183_G2 mul_by_q() const;
    bn183_G2 mul_by_cofactor() const;

    bool is_well_formed() const;

    static bn183_G2 zero();
    static bn183_G2 one();
    static bn183_G2 random_element();

    static std::size_t size_in_bits() { return twist_field::ceil_size_in_bits() + 1; }
    static bigint<base_field::num_limbs> field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const bn183_G2 &g);
    friend std::istream& operator>>(std::istream &in, bn183_G2 &g);

    static void batch_to_special_all_non_zeros(std::vector<bn183_G2> &vec);
};

template<mp_size_t m>
bn183_G2 operator*(const bigint<m> &lhs, const bn183_G2 &rhs)
{
    return scalar_mul<bn183_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
bn183_G2 operator*(const Fp_model<m,modulus_p> &lhs, const bn183_G2 &rhs)
{
    return scalar_mul<bn183_G2, m>(rhs, lhs.as_bigint());
}


} // namespace libff
#endif // BN183_G2_HPP_
