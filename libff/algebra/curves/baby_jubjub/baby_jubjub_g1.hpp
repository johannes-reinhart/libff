#ifndef BABY_JUBJUB_G1_HPP_
#define BABY_JUBJUB_G1_HPP_
#include <vector>

#include <libff/algebra/curves/baby_jubjub/baby_jubjub_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class baby_jubjub_G1;
std::ostream& operator<<(std::ostream &, const baby_jubjub_G1&);
std::istream& operator>>(std::istream &, baby_jubjub_G1&);

class baby_jubjub_G1 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static baby_jubjub_G1 G1_zero;
    static baby_jubjub_G1 G1_one;
    static baby_jubjub_Fq coeff_a;
    static baby_jubjub_Fq coeff_d;
    static bool initialized;

    baby_jubjub_Fq X, Y, Z;
    baby_jubjub_G1();
private:
    baby_jubjub_G1(const baby_jubjub_Fq& X, const baby_jubjub_Fq& Y, const baby_jubjub_Fq& Z) : X(X), Y(Y), Z(Z) {};

public:
    typedef baby_jubjub_Fq base_field;
    typedef baby_jubjub_Fr scalar_field;
    // using inverted coordinates
    baby_jubjub_G1(const baby_jubjub_Fq& X, const baby_jubjub_Fq& Y) : X(Y), Y(X), Z(X*Y) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const baby_jubjub_G1 &other) const;
    bool operator!=(const baby_jubjub_G1 &other) const;

    baby_jubjub_G1 operator+(const baby_jubjub_G1 &other) const;
    baby_jubjub_G1 operator-() const;
    baby_jubjub_G1 operator-(const baby_jubjub_G1 &other) const;

    baby_jubjub_G1 add(const baby_jubjub_G1 &other) const;
    baby_jubjub_G1 mixed_add(const baby_jubjub_G1 &other) const;
    baby_jubjub_G1 dbl() const;

    static baby_jubjub_G1 from_y(const baby_jubjub_Fq &Y);

    bool is_well_formed() const;

    static baby_jubjub_G1 zero();
    static baby_jubjub_G1 one();
    static baby_jubjub_G1 random_element();

    static std::size_t size_in_bits() { return baby_jubjub_Fq::ceil_size_in_bits() + 1; }
    static bigint<base_field::num_limbs> field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const baby_jubjub_G1 &g);
    friend std::istream& operator>>(std::istream &in, baby_jubjub_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<baby_jubjub_G1> &vec);
};

template<mp_size_t m>
baby_jubjub_G1 operator*(const bigint<m> &lhs, const baby_jubjub_G1 &rhs)
{
    return scalar_mul<baby_jubjub_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
baby_jubjub_G1 operator*(const Fp_model<m,modulus_p> &lhs, const baby_jubjub_G1 &rhs)
{
    return scalar_mul<baby_jubjub_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<baby_jubjub_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<baby_jubjub_G1> &v);

} // libff
#endif // BABY_JUBJUB_G1_HPP_
