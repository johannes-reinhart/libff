#include <libff/algebra/curves/bn254/bn254_g2.hpp>

namespace libff {

using std::size_t;

#ifdef PROFILE_OP_COUNTS
long long bn254_G2::add_cnt = 0;
long long bn254_G2::dbl_cnt = 0;
#endif

std::vector<size_t> bn254_G2::wnaf_window_table;
std::vector<size_t> bn254_G2::fixed_base_exp_window_table;
bn254_G2 bn254_G2::G2_zero = {};
bn254_G2 bn254_G2::G2_one = {};
bool bn254_G2::initialized = false;
bigint<bn254_G2::h_limbs> bn254_G2::h;

bn254_G2::bn254_G2()
{
    if (initialized)
    {
        this->X = G2_zero.X;
        this->Y = G2_zero.Y;
        this->Z = G2_zero.Z;
    }
}

bn254_Fq2 bn254_G2::mul_by_b(const bn254_Fq2 &elt)
{
    return bn254_Fq2(bn254_twist_mul_by_b_c0 * elt.c0, bn254_twist_mul_by_b_c1 * elt.c1);
}

void bn254_G2::print() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        bn254_G2 copy(*this);
        copy.to_affine_coordinates();
        gmp_printf("(%Nd*z + %Nd , %Nd*z + %Nd)\n",
                   copy.X.c1.as_bigint().data, bn254_Fq::num_limbs,
                   copy.X.c0.as_bigint().data, bn254_Fq::num_limbs,
                   copy.Y.c1.as_bigint().data, bn254_Fq::num_limbs,
                   copy.Y.c0.as_bigint().data, bn254_Fq::num_limbs);
    }
}

void bn254_G2::print_coordinates() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        gmp_printf("(%Nd*z + %Nd : %Nd*z + %Nd : %Nd*z + %Nd)\n",
                   this->X.c1.as_bigint().data, bn254_Fq::num_limbs,
                   this->X.c0.as_bigint().data, bn254_Fq::num_limbs,
                   this->Y.c1.as_bigint().data, bn254_Fq::num_limbs,
                   this->Y.c0.as_bigint().data, bn254_Fq::num_limbs,
                   this->Z.c1.as_bigint().data, bn254_Fq::num_limbs,
                   this->Z.c0.as_bigint().data, bn254_Fq::num_limbs);
    }
}

void bn254_G2::to_affine_coordinates()
{
    if (this->is_zero())
    {
        this->X = bn254_Fq2::zero();
        this->Y = bn254_Fq2::one();
        this->Z = bn254_Fq2::zero();
    }
    else
    {
        bn254_Fq2 Z_inv = Z.inverse();
        bn254_Fq2 Z2_inv = Z_inv.squared();
        bn254_Fq2 Z3_inv = Z2_inv * Z_inv;
        this->X = this->X * Z2_inv;
        this->Y = this->Y * Z3_inv;
        this->Z = bn254_Fq2::one();
    }
}

void bn254_G2::to_special()
{
    this->to_affine_coordinates();
}

bool bn254_G2::is_special() const
{
    return (this->is_zero() || this->Z == bn254_Fq2::one());
}

bool bn254_G2::is_zero() const
{
    return (this->Z.is_zero());
}

bool bn254_G2::operator==(const bn254_G2 &other) const
{
    if (this->is_zero())
    {
        return other.is_zero();
    }

    if (other.is_zero())
    {
        return false;
    }

    /* now neither is O */

    // using Jacobian coordinates so:
    // (X1:Y1:Z1) = (X2:Y2:Z2)
    // iff
    // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    // iff
    // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

    bn254_Fq2 Z1_squared = (this->Z).squared();
    bn254_Fq2 Z2_squared = (other.Z).squared();

    if ((this->X * Z2_squared) != (other.X * Z1_squared))
    {
        return false;
    }

    bn254_Fq2 Z1_cubed = (this->Z) * Z1_squared;
    bn254_Fq2 Z2_cubed = (other.Z) * Z2_squared;

    return !((this->Y * Z2_cubed) != (other.Y * Z1_cubed));
}

bool bn254_G2::operator!=(const bn254_G2& other) const
{
    return !(operator==(other));
}

bn254_G2 bn254_G2::operator+(const bn254_G2 &other) const
{
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // using Jacobian coordinates according to
    // https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
    // Note: (X1:Y1:Z1) = (X2:Y2:Z2)
    // iff
    // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    // iff
    // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

    bn254_Fq2 Z1Z1 = (this->Z).squared();
    bn254_Fq2 Z2Z2 = (other.Z).squared();

    bn254_Fq2 U1 = this->X * Z2Z2;
    bn254_Fq2 U2 = other.X * Z1Z1;

    bn254_Fq2 Z1_cubed = (this->Z) * Z1Z1;
    bn254_Fq2 Z2_cubed = (other.Z) * Z2Z2;

    bn254_Fq2 S1 = (this->Y) * Z2_cubed;      // S1 = Y1 * Z2 * Z2Z2
    bn254_Fq2 S2 = (other.Y) * Z1_cubed;      // S2 = Y2 * Z1 * Z1Z1

    // check for doubling case
    if (U1 == U2 && S1 == S2)
    {
        // dbl case; nothing of above can be reused
        return this->dbl();
    }

#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif

    // rest of add case
    bn254_Fq2 H = U2 - U1;                            // H = U2-U1
    bn254_Fq2 S2_minus_S1 = S2-S1;
    bn254_Fq2 I = (H+H).squared();                    // I = (2 * H)^2
    bn254_Fq2 J = H * I;                              // J = H * I
    bn254_Fq2 r = S2_minus_S1 + S2_minus_S1;          // r = 2 * (S2-S1)
    bn254_Fq2 V = U1 * I;                             // V = U1 * I
    bn254_Fq2 X3 = r.squared() - J - (V+V);           // X3 = r^2 - J - 2 * V
    bn254_Fq2 S1_J = S1 * J;
    bn254_Fq2 Y3 = r * (V-X3) - (S1_J+S1_J);          // Y3 = r * (V-X3)-2 S1 J
    bn254_Fq2 Z3 = ((this->Z+other.Z).squared()-Z1Z1-Z2Z2) * H; // Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2) * H

    return bn254_G2(X3, Y3, Z3);
}

bn254_G2 bn254_G2::operator-() const
{
    return bn254_G2(this->X, -(this->Y), this->Z);
}


bn254_G2 bn254_G2::operator-(const bn254_G2 &other) const
{
    return (*this) + (-other);
}

bn254_G2 bn254_G2::add(const bn254_G2 &other) const
{
       return (*this) + other;
}

bn254_G2 bn254_G2::mixed_add(const bn254_G2 &other) const
{
#ifdef DEBUG
    assert(other.is_special());
#endif

    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // using Jacobian coordinates according to
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
    // Note: (X1:Y1:Z1) = (X2:Y2:Z2)
    // iff
    // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    // iff
    // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3
    // we know that Z2 = 1

    const bn254_Fq2 Z1Z1 = (this->Z).squared();

    const bn254_Fq2 &U1 = this->X;
    const bn254_Fq2 U2 = other.X * Z1Z1;

    const bn254_Fq2 Z1_cubed = (this->Z) * Z1Z1;

    const bn254_Fq2 &S1 = (this->Y);                // S1 = Y1 * Z2 * Z2Z2
    const bn254_Fq2 S2 = (other.Y) * Z1_cubed;      // S2 = Y2 * Z1 * Z1Z1

    // check for doubling case
    if (U1 == U2 && S1 == S2)
    {
        // dbl case; nothing of above can be reused
        return this->dbl();
    }

#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif

    bn254_Fq2 H = U2-(this->X);                         // H = U2-X1
    bn254_Fq2 HH = H.squared() ;                        // HH = H^2
    bn254_Fq2 I = HH+HH;                                // I = 4*HH
    I = I + I;
    bn254_Fq2 J = H*I;                                  // J = H*I
    bn254_Fq2 r = S2-(this->Y);                         // r = 2*(S2-Y1)
    r = r + r;
    bn254_Fq2 V = (this->X) * I ;                       // V = X1*I
    bn254_Fq2 X3 = r.squared()-J-V-V;                   // X3 = r^2-J-2*V
    bn254_Fq2 Y3 = (this->Y)*J;                         // Y3 = r*(V-X3)-2*Y1*J
    Y3 = r*(V-X3) - Y3 - Y3;
    bn254_Fq2 Z3 = ((this->Z)+H).squared() - Z1Z1 - HH; // Z3 = (Z1+H)^2-Z1Z1-HH

    return bn254_G2(X3, Y3, Z3);
}

bn254_G2 bn254_G2::dbl() const
{
#ifdef PROFILE_OP_COUNTS
    this->dbl_cnt++;
#endif
    // handle point at infinity
    if (this->is_zero())
    {
        return (*this);
    }

    // NOTE: does not handle O and pts of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // using Jacobian coordinates according to
    // https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
    bn254_Fq2 A = (this->X).squared();         // A = X1^2
    bn254_Fq2 B = (this->Y).squared();        // B = Y1^2
    bn254_Fq2 C = B.squared();                // C = B^2
    bn254_Fq2 D = (this->X + B).squared() - A - C;
    D = D+D;                        // D = 2 * ((X1 + B)^2 - A - C)
    bn254_Fq2 E = A + A + A;                  // E = 3 * A
    bn254_Fq2 F = E.squared();                // F = E^2
    bn254_Fq2 X3 = F - (D+D);                 // X3 = F - 2 D
    bn254_Fq2 eightC = C+C;
    eightC = eightC + eightC;
    eightC = eightC + eightC;
    bn254_Fq2 Y3 = E * (D - X3) - eightC;     // Y3 = E * (D - X3) - 8 * C
    bn254_Fq2 Y1Z1 = (this->Y)*(this->Z);
    bn254_Fq2 Z3 = Y1Z1 + Y1Z1;               // Z3 = 2 * Y1 * Z1

    return bn254_G2(X3, Y3, Z3);
}

bn254_G2 bn254_G2::mul_by_q() const
{
    return bn254_G2(bn254_twist_mul_by_q_X * (this->X).Frobenius_map(1),
                      bn254_twist_mul_by_q_Y * (this->Y).Frobenius_map(1),
                      (this->Z).Frobenius_map(1));
}

bn254_G2 bn254_G2::mul_by_cofactor() const
{
    return bn254_G2::h * (*this);
}

bool bn254_G2::is_well_formed() const
{
    if (this->is_zero())
    {
        return true;
    }
    /*
      y^2 = x^3 + b

      We are using Jacobian coordinates, so equation we need to check is actually

      (y/z^3)^2 = (x/z^2)^3 + b
      y^2 / z^6 = x^3 / z^6 + b
      y^2 = x^3 + b z^6
    */
    bn254_Fq2 X2 = this->X.squared();
    bn254_Fq2 Y2 = this->Y.squared();
    bn254_Fq2 Z2 = this->Z.squared();

    bn254_Fq2 X3 = this->X * X2;
    bn254_Fq2 Z3 = this->Z * Z2;
    bn254_Fq2 Z6 = Z3.squared();

    return (Y2 == X3 + bn254_twist_coeff_b * Z6);
}

bn254_G2 bn254_G2::zero()
{
    return G2_zero;
}

bn254_G2 bn254_G2::one()
{
    return G2_one;
}

bn254_G2 bn254_G2::random_element()
{
    return (bn254_Fr::random_element().as_bigint()) * G2_one;
}

std::ostream& operator<<(std::ostream &out, const bn254_G2 &g)
{
    bn254_G2 copy(g);
    copy.to_affine_coordinates();
    out << (copy.is_zero() ? 1 : 0) << OUTPUT_SEPARATOR;
#ifdef NO_PT_COMPRESSION
    out << copy.X << OUTPUT_SEPARATOR << copy.Y;
#else
    /* storing LSB of Y */
    out << copy.X << OUTPUT_SEPARATOR << (copy.Y.c0.as_bigint().data[0] & 1);
#endif

    return out;
}

std::istream& operator>>(std::istream &in, bn254_G2 &g)
{
    char is_zero;
    bn254_Fq2 tX, tY;

#ifdef NO_PT_COMPRESSION
    in >> is_zero >> tX >> tY;
    is_zero -= '0';
#else
    in.read((char*)&is_zero, 1); // this reads is_zero;
    is_zero -= '0';
    consume_OUTPUT_SEPARATOR(in);

    unsigned char Y_lsb;
    in >> tX;
    consume_OUTPUT_SEPARATOR(in);
    in.read((char*)&Y_lsb, 1);
    Y_lsb -= '0';

    // y = +/- sqrt(x^3 + b)
    if (is_zero == 0)
    {
        bn254_Fq2 tX2 = tX.squared();
        bn254_Fq2 tY2 = tX2 * tX + bn254_twist_coeff_b;
        tY = tY2.sqrt();

        if ((tY.c0.as_bigint().data[0] & 1) != Y_lsb)
        {
            tY = -tY;
        }
    }
#endif
    // using projective coordinates
    if (is_zero == 0)
    {
        g.X = tX;
        g.Y = tY;
        g.Z = bn254_Fq2::one();
    }
    else
    {
        g = bn254_G2::zero();
    }

    return in;
}

void bn254_G2::batch_to_special_all_non_zeros(std::vector<bn254_G2> &vec)
{
    std::vector<bn254_Fq2> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el: vec)
    {
        Z_vec.emplace_back(el.Z);
    }
    batch_invert<bn254_Fq2>(Z_vec);

    const bn254_Fq2 one = bn254_Fq2::one();

    for (size_t i = 0; i < vec.size(); ++i)
    {
        bn254_Fq2 Z2 = Z_vec[i].squared();
        bn254_Fq2 Z3 = Z_vec[i] * Z2;

        vec[i].X = vec[i].X * Z2;
        vec[i].Y = vec[i].Y * Z3;
        vec[i].Z = one;
    }
}

} // namespace libff
