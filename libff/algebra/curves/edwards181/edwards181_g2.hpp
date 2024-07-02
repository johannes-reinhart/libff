#ifndef EDWARDS181_G2_HPP_
#define EDWARDS181_G2_HPP_
#include <vector>

#include <libff/algebra/curves/edwards181/edwards181_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class edwards181_G2;
std::ostream& operator<<(std::ostream &, const edwards181_G2&);
std::istream& operator>>(std::istream &, edwards181_G2&);

class edwards181_G2 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static edwards181_G2 G2_zero;
    static edwards181_G2 G2_one;
    static bool initialized;

    edwards181_Fq3 X, Y, Z;
    edwards181_G2();
private:
    edwards181_G2(const edwards181_Fq3& X, const edwards181_Fq3& Y, const edwards181_Fq3& Z) : X(X), Y(Y), Z(Z) {};
public:
    static edwards181_Fq3 mul_by_twist(const edwards181_Fq3 &elt);
    static edwards181_Fq3 mul_by_a(const edwards181_Fq3 &elt);
    static edwards181_Fq3 mul_by_d(const edwards181_Fq3 &elt);
    typedef edwards181_Fq base_field;
    typedef edwards181_Fq3 twist_field;
    typedef edwards181_Fr scalar_field;

    // using inverted coordinates
    edwards181_G2(const edwards181_Fq3& X, const edwards181_Fq3& Y) : X(Y), Y(X), Z(X*Y) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const edwards181_G2 &other) const;
    bool operator!=(const edwards181_G2 &other) const;

    edwards181_G2 operator+(const edwards181_G2 &other) const;
    edwards181_G2 operator-() const;
    edwards181_G2 operator-(const edwards181_G2 &other) const;

    edwards181_G2 add(const edwards181_G2 &other) const;
    edwards181_G2 mixed_add(const edwards181_G2 &other) const;
    edwards181_G2 dbl() const;
    edwards181_G2 mul_by_q() const;

    bool is_well_formed() const;

    static edwards181_G2 zero();
    static edwards181_G2 one();
    static edwards181_G2 random_element();

    static size_t size_in_bits() { return twist_field::ceil_size_in_bits() + 1; }
    static bigint<base_field::num_limbs> field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const edwards181_G2 &g);
    friend std::istream& operator>>(std::istream &in, edwards181_G2 &g);

    static void batch_to_special_all_non_zeros(std::vector<edwards181_G2> &vec);
};

template<mp_size_t m>
edwards181_G2 operator*(const bigint<m> &lhs, const edwards181_G2 &rhs)
{
    return scalar_mul<edwards181_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
edwards181_G2 operator*(const Fp_model<m,modulus_p> &lhs, const edwards181_G2 &rhs)
{
    return scalar_mul<edwards181_G2, m>(rhs, lhs.as_bigint());
}


} // libff
#endif // EDWARDS181_G2_HPP_
