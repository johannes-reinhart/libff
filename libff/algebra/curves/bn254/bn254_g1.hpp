#ifndef BN254_G1_HPP_
#define BN254_G1_HPP_
#include <vector>

#include <libff/algebra/curves/bn254/bn254_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class bn254_G1;
std::ostream& operator<<(std::ostream &, const bn254_G1&);
std::istream& operator>>(std::istream &, bn254_G1&);

class bn254_G1 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<std::size_t> wnaf_window_table;
    static std::vector<std::size_t> fixed_base_exp_window_table;
    static bn254_G1 G1_zero;
    static bn254_G1 G1_one;
    static bool initialized;

    typedef bn254_Fq base_field;
    typedef bn254_Fr scalar_field;

    // Cofactor
    static const mp_size_t h_bitcount = 1;
    static const mp_size_t h_limbs = (h_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
    static bigint<h_limbs> h;

    bn254_Fq X, Y, Z;

    // using Jacobian coordinates
    bn254_G1();
    bn254_G1(const bn254_Fq& X, const bn254_Fq& Y) : X(X), Y(Y), Z(1) {};
    bn254_G1(const bn254_Fq& X, const bn254_Fq& Y, const bn254_Fq& Z) : X(X), Y(Y), Z(Z) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const bn254_G1 &other) const;
    bool operator!=(const bn254_G1 &other) const;

    bn254_G1 operator+(const bn254_G1 &other) const;
    bn254_G1 operator-() const;
    bn254_G1 operator-(const bn254_G1 &other) const;

    bn254_G1 add(const bn254_G1 &other) const;
    bn254_G1 mixed_add(const bn254_G1 &other) const;
    bn254_G1 dbl() const;
    bn254_G1 mul_by_cofactor() const;

    bool is_well_formed() const;

    static bn254_G1 zero();
    static bn254_G1 one();
    static bn254_G1 random_element();

    static std::size_t size_in_bits() { return base_field::ceil_size_in_bits() + 1; }
    static bigint<base_field::num_limbs> field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const bn254_G1 &g);
    friend std::istream& operator>>(std::istream &in, bn254_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<bn254_G1> &vec);
};

template<mp_size_t m>
bn254_G1 operator*(const bigint<m> &lhs, const bn254_G1 &rhs)
{
    return scalar_mul<bn254_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
bn254_G1 operator*(const Fp_model<m,modulus_p> &lhs, const bn254_G1 &rhs)
{
    return scalar_mul<bn254_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<bn254_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<bn254_G1> &v);

} // namespace libff
#endif // BN254_G1_HPP_
