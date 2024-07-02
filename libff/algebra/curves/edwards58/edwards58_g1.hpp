#ifndef EDWARDS58_G1_HPP_
#define EDWARDS58_G1_HPP_
#include <vector>

#include <libff/algebra/curves/edwards58/edwards58_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class edwards58_G1;
std::ostream& operator<<(std::ostream &, const edwards58_G1&);
std::istream& operator>>(std::istream &, edwards58_G1&);

class edwards58_G1 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static edwards58_G1 G1_zero;
    static edwards58_G1 G1_one;
    static edwards58_Fq coeff_a;
    static edwards58_Fq coeff_d;
    static bool initialized;

    edwards58_Fq X, Y, Z;
    edwards58_G1();
private:
    edwards58_G1(const edwards58_Fq& X, const edwards58_Fq& Y, const edwards58_Fq& Z) : X(X), Y(Y), Z(Z) {};

public:
    typedef edwards58_Fq base_field;
    typedef edwards58_Fr scalar_field;
    // using inverted coordinates
    edwards58_G1(const edwards58_Fq& X, const edwards58_Fq& Y) : X(Y), Y(X), Z(X*Y) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const edwards58_G1 &other) const;
    bool operator!=(const edwards58_G1 &other) const;

    edwards58_G1 operator+(const edwards58_G1 &other) const;
    edwards58_G1 operator-() const;
    edwards58_G1 operator-(const edwards58_G1 &other) const;

    edwards58_G1 add(const edwards58_G1 &other) const;
    edwards58_G1 mixed_add(const edwards58_G1 &other) const;
    edwards58_G1 dbl() const;

    static edwards58_G1 from_y(const edwards58_Fq &Y);

    bool is_well_formed() const;

    static edwards58_G1 zero();
    static edwards58_G1 one();
    static edwards58_G1 random_element();

    static std::size_t size_in_bits() { return edwards58_Fq::ceil_size_in_bits() + 1; }
    static bigint<base_field::num_limbs> field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const edwards58_G1 &g);
    friend std::istream& operator>>(std::istream &in, edwards58_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<edwards58_G1> &vec);
};

template<mp_size_t m>
edwards58_G1 operator*(const bigint<m> &lhs, const edwards58_G1 &rhs)
{
    return scalar_mul<edwards58_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
edwards58_G1 operator*(const Fp_model<m,modulus_p> &lhs, const edwards58_G1 &rhs)
{
    return scalar_mul<edwards58_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<edwards58_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<edwards58_G1> &v);

} // libff
#endif // EDWARDS58_G1_HPP_
