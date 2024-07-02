#ifndef JUBJUB_BN183_G1_HPP_
#define JUBJUB_BN183_G1_HPP_
#include <vector>

#include <libff/algebra/curves/jubjub_bn183/jubjub_bn183_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class jubjub_bn183_G1;
std::ostream& operator<<(std::ostream &, const jubjub_bn183_G1&);
std::istream& operator>>(std::istream &, jubjub_bn183_G1&);

class jubjub_bn183_G1 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static jubjub_bn183_G1 G1_zero;
    static jubjub_bn183_G1 G1_one;
    static jubjub_bn183_Fq coeff_a;
    static jubjub_bn183_Fq coeff_d;
    static bool initialized;

    jubjub_bn183_Fq X, Y, Z;
    jubjub_bn183_G1();
private:
    jubjub_bn183_G1(const jubjub_bn183_Fq& X, const jubjub_bn183_Fq& Y, const jubjub_bn183_Fq& Z) : X(X), Y(Y), Z(Z) {};

public:
    typedef jubjub_bn183_Fq base_field;
    typedef jubjub_bn183_Fr scalar_field;
    // using inverted coordinates
    jubjub_bn183_G1(const jubjub_bn183_Fq& X, const jubjub_bn183_Fq& Y) : X(Y), Y(X), Z(X*Y) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const jubjub_bn183_G1 &other) const;
    bool operator!=(const jubjub_bn183_G1 &other) const;

    jubjub_bn183_G1 operator+(const jubjub_bn183_G1 &other) const;
    jubjub_bn183_G1 operator-() const;
    jubjub_bn183_G1 operator-(const jubjub_bn183_G1 &other) const;

    jubjub_bn183_G1 add(const jubjub_bn183_G1 &other) const;
    jubjub_bn183_G1 mixed_add(const jubjub_bn183_G1 &other) const;
    jubjub_bn183_G1 dbl() const;

    static jubjub_bn183_G1 from_y(const jubjub_bn183_Fq &Y);

    bool is_well_formed() const;

    static jubjub_bn183_G1 zero();
    static jubjub_bn183_G1 one();
    static jubjub_bn183_G1 random_element();

    static std::size_t size_in_bits() { return jubjub_bn183_Fq::ceil_size_in_bits() + 1; }
    static bigint<base_field::num_limbs> field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const jubjub_bn183_G1 &g);
    friend std::istream& operator>>(std::istream &in, jubjub_bn183_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<jubjub_bn183_G1> &vec);
};

template<mp_size_t m>
jubjub_bn183_G1 operator*(const bigint<m> &lhs, const jubjub_bn183_G1 &rhs)
{
    return scalar_mul<jubjub_bn183_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
jubjub_bn183_G1 operator*(const Fp_model<m,modulus_p> &lhs, const jubjub_bn183_G1 &rhs)
{
    return scalar_mul<jubjub_bn183_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<jubjub_bn183_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<jubjub_bn183_G1> &v);

} // libff
#endif // JUBJUB_BN183_G1_HPP_
