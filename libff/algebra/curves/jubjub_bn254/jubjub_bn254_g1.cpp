#include <libff/algebra/curves/jubjub_bn254/jubjub_bn254_g1.hpp>

namespace libff {

#ifdef PROFILE_OP_COUNTS
long long jubjub_bn254_G1::add_cnt = 0;
long long jubjub_bn254_G1::dbl_cnt = 0;
#endif

std::vector<size_t> jubjub_bn254_G1::wnaf_window_table;
std::vector<size_t> jubjub_bn254_G1::fixed_base_exp_window_table;
jubjub_bn254_G1 jubjub_bn254_G1::G1_zero;
jubjub_bn254_G1 jubjub_bn254_G1::G1_one;
jubjub_bn254_Fq jubjub_bn254_G1::coeff_a;
jubjub_bn254_Fq jubjub_bn254_G1::coeff_d;
bool jubjub_bn254_G1::initialized = false;

jubjub_bn254_G1::jubjub_bn254_G1()
{
    if(initialized)
    {
        this->X = G1_zero.X;
        this->Y = G1_zero.Y;
        this->Z = G1_zero.Z;
    }
}

void jubjub_bn254_G1::print() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        jubjub_bn254_G1 copy(*this);
        copy.to_affine_coordinates();
        gmp_printf("(%Nd , %Nd)\n",
                   copy.X.as_bigint().data, jubjub_bn254_Fq::num_limbs,
                   copy.Y.as_bigint().data, jubjub_bn254_Fq::num_limbs);
    }
}

void jubjub_bn254_G1::print_coordinates() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        gmp_printf("(%Nd : %Nd : %Nd)\n",
                   this->X.as_bigint().data, jubjub_bn254_Fq::num_limbs,
                   this->Y.as_bigint().data, jubjub_bn254_Fq::num_limbs,
                   this->Z.as_bigint().data, jubjub_bn254_Fq::num_limbs);
    }
}

void jubjub_bn254_G1::to_affine_coordinates()
{
    if (this->is_zero())
    {
        this->X = jubjub_bn254_Fq::zero();
        this->Y = jubjub_bn254_Fq::one();
        this->Z = jubjub_bn254_Fq::one();
    }
    else
    {
        // go from inverted coordinates to projective coordinates
        jubjub_bn254_Fq tX = this->Y * this->Z;
        jubjub_bn254_Fq tY = this->X * this->Z;
        jubjub_bn254_Fq tZ = this->X * this->Y;
        // go from projective coordinates to affine coordinates
        jubjub_bn254_Fq tZ_inv = tZ.inverse();
        this->X = tX * tZ_inv;
        this->Y = tY * tZ_inv;
        this->Z = jubjub_bn254_Fq::one();
    }
}


jubjub_bn254_G1 jubjub_bn254_G1::from_y(const jubjub_bn254_Fq &Y){
    for (int i = 0; i < 1000; i++){
        const jubjub_bn254_Fq y = Y+i;
        const jubjub_bn254_Fq yy = y.squared();
        const jubjub_bn254_Fq lhs = (yy - 1);
        const jubjub_bn254_Fq rhs = (coeff_d*yy - coeff_a);
        const jubjub_bn254_Fq xx = lhs * rhs.inverse();

        if (jacobi(xx) == -1){
            continue;
        }

        jubjub_bn254_Fq x = xx.sqrt();
        // return point with negative x coordinate (Ethsnarks defines negative other way than I think is intuitive? TODO)
        if (is_negative(x)){
            return jubjub_bn254_G1(x, y);
        }else{
            return jubjub_bn254_G1(-x, y);
        }
    }

    throw std::runtime_error("Could not find a point on curve");
}


void jubjub_bn254_G1::to_special()
{
    if (this->Z.is_zero())
    {
        return;
    }

#ifdef DEBUG
    const jubjub_bn254_G1 copy(*this);
#endif

    jubjub_bn254_Fq Z_inv = this->Z.inverse();
    this->X = this->X * Z_inv;
    this->Y = this->Y * Z_inv;
    this->Z = jubjub_bn254_Fq::one();

#ifdef DEBUG
    assert((*this) == copy);
#endif
}

bool jubjub_bn254_G1::is_special() const
{
    return (this->is_zero() || this->Z == jubjub_bn254_Fq::one());
}

bool jubjub_bn254_G1::is_zero() const
{
    return (this->Y.is_zero() && this->Z.is_zero());
}

bool jubjub_bn254_G1::operator==(const jubjub_bn254_G1 &other) const
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

    // X1/Z1 = X2/Z2 <=> X1*Z2 = X2*Z1
    if ((this->X * other.Z) != (other.X * this->Z))
    {
        return false;
    }

    // Y1/Z1 = Y2/Z2 <=> Y1*Z2 = Y2*Z1
    if ((this->Y * other.Z) != (other.Y * this->Z))
    {
        return false;
    }

    return true;
}

bool jubjub_bn254_G1::operator!=(const jubjub_bn254_G1& other) const
{
    return !(operator==(other));
}

jubjub_bn254_G1 jubjub_bn254_G1::operator+(const jubjub_bn254_G1 &other) const
{
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return (*this);
    }

    return this->add(other);
}

jubjub_bn254_G1 jubjub_bn254_G1::operator-() const
{
    return jubjub_bn254_G1(-(this->X), this->Y, this->Z);
}


jubjub_bn254_G1 jubjub_bn254_G1::operator-(const jubjub_bn254_G1 &other) const
{
    return (*this) + (-other);
}

jubjub_bn254_G1 jubjub_bn254_G1::add(const jubjub_bn254_G1 &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-twisted-inverted.html#addition-add-2008-bbjlp

    jubjub_bn254_Fq A = (this->Z) * (other.Z);                   // A = Z1*Z2
    jubjub_bn254_Fq B = coeff_d * A.squared();    // B = d*A^2
    jubjub_bn254_Fq C = (this->X) * (other.X);                   // C = X1*X2
    jubjub_bn254_Fq D = (this->Y) * (other.Y);                   // D = Y1*Y2
    jubjub_bn254_Fq E = C * D;                                   // E = C*D
    jubjub_bn254_Fq H = C - coeff_a*D;            // H = C-a*D
    jubjub_bn254_Fq I = (this->X+this->Y)*(other.X+other.Y)-C-D; // I = (X1+Y1)*(X2+Y2)-C-D
    jubjub_bn254_Fq X3 = (E+B)*H;                                // X3 = c*(E+B)*H
    jubjub_bn254_Fq Y3 = (E-B)*I;                                // Y3 = c*(E-B)*I
    jubjub_bn254_Fq Z3 = A*H*I;                                  // Z3 = A*H*I

    return jubjub_bn254_G1(X3, Y3, Z3);
}

jubjub_bn254_G1 jubjub_bn254_G1::mixed_add(const jubjub_bn254_G1 &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
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

#ifdef DEBUG
    assert(other.is_special());
#endif

    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-twisted-inverted.html#addition-madd-2008-bbjlp

    jubjub_bn254_Fq A = this->Z;                                 // A = Z1
    jubjub_bn254_Fq B = coeff_d * A.squared();    // B = d*A^2
    jubjub_bn254_Fq C = (this->X) * (other.X);                   // C = X1*X2
    jubjub_bn254_Fq D = (this->Y) * (other.Y);                   // D = Y1*Y2
    jubjub_bn254_Fq E = C * D;                                   // E = C*D
    jubjub_bn254_Fq H = C - coeff_a*D;                  // H = C-a*D
    jubjub_bn254_Fq I = (this->X+this->Y)*(other.X+other.Y)-C-D; // I = (X1+Y1)*(X2+Y2)-C-D
    jubjub_bn254_Fq X3 = (E+B)*H;                                // X3 = c*(E+B)*H
    jubjub_bn254_Fq Y3 = (E-B)*I;                                // Y3 = c*(E-B)*I
    jubjub_bn254_Fq Z3 = A*H*I;                                  // Z3 = A*H*I

    return jubjub_bn254_G1(X3, Y3, Z3);
}

jubjub_bn254_G1 jubjub_bn254_G1::dbl() const
{
#ifdef PROFILE_OP_COUNTS
    this->dbl_cnt++;
#endif
    if (this->is_zero())
    {
        return (*this);
    }
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-twisted-inverted.html#doubling-dbl-2008-bbjlp

    jubjub_bn254_Fq A = (this->X).squared();                      // A = X1^2
    jubjub_bn254_Fq B = (this->Y).squared();                      // B = Y1^2
    jubjub_bn254_Fq U = coeff_a*B;                    // U = a*B
    jubjub_bn254_Fq C = A+U;                                      // C = A+U
    jubjub_bn254_Fq D = A-U;                                      // D = A-U
    jubjub_bn254_Fq E = (this->X+this->Y).squared()-A-B;          // E = (X1+Y1)^2-A-B
    jubjub_bn254_Fq X3 = C*D;                                     // X3 = C*D
    jubjub_bn254_Fq dZZ = coeff_d * this->Z.squared();
    jubjub_bn254_Fq Y3 = E*(C-dZZ-dZZ);                           // Y3 = E*(C-2*d*Z1^2)
    jubjub_bn254_Fq Z3 = D*E;                                     // Z3 = D*E

    return jubjub_bn254_G1(X3, Y3, Z3);
}

bool jubjub_bn254_G1::is_well_formed() const
{
     /* Note that point at infinity is the only special case we must check as
       inverted representation does no cover points (0, +-c) and (+-c, 0). */
    if (this->is_zero())
    {
        return true;
    }
    /*
        a x^2 + y^2 = 1 + d x^2 y^2

        We are using inverted, so equation we need to check is actually

        a (z/x)^2 + (z/y)^2 = 1 + d z^4 / (x^2 * y^2)
        z^2 (a y^2 + x^2 - dz^2) = x^2 y^2
    */
    jubjub_bn254_Fq X2 = this->X.squared();
    jubjub_bn254_Fq Y2 = this->Y.squared();
    jubjub_bn254_Fq Z2 = this->Z.squared();


    return (Z2 * (coeff_a * Y2 + X2 - coeff_d * Z2) == X2 * Y2);
}

jubjub_bn254_G1 jubjub_bn254_G1::zero()
{
    return G1_zero;
}

jubjub_bn254_G1 jubjub_bn254_G1::one()
{
    return G1_one;
}

jubjub_bn254_G1 jubjub_bn254_G1::random_element()
{
    return jubjub_bn254_Fr::random_element().as_bigint() * G1_one;
}

std::ostream& operator<<(std::ostream &out, const jubjub_bn254_G1 &g)
{
    jubjub_bn254_G1 copy(g);
    copy.to_affine_coordinates();
#ifdef NO_PT_COMPRESSION
    out << copy.X << OUTPUT_SEPARATOR << copy.Y;
#else
    /* storing LSB of Y */
    out << copy.X << OUTPUT_SEPARATOR << (copy.Y.as_bigint().data[0] & 1);
#endif

    return out;
}

std::istream& operator>>(std::istream &in, jubjub_bn254_G1 &g)
{
    jubjub_bn254_Fq tX, tY;

#ifdef NO_PT_COMPRESSION
    in >> tX;
    consume_OUTPUT_SEPARATOR(in);
    in >> tY;
#else
    /*
      a x^2 + y^2 = 1 + d x^2 y^2
      y = sqrt((1-ax^2)/(1-dx^2))
    */
    unsigned char Y_lsb;
    in >> tX;

    consume_OUTPUT_SEPARATOR(in);
    in.read((char*)&Y_lsb, 1);
    Y_lsb -= '0';

    jubjub_bn254_Fq tX2 = tX.squared();
    jubjub_bn254_Fq tY2 = (jubjub_bn254_Fq::one() - jubjub_bn254_G1::coeff_a * tX2) *
        (jubjub_bn254_Fq::one() - jubjub_bn254_G1::coeff_d * tX2).inverse();
    tY = tY2.sqrt();

    if ((tY.as_bigint().data[0] & 1) != Y_lsb)
    {
        tY = -tY;
    }
#endif

    // using inverted coordinates
    g.X = tY;
    g.Y = tX;
    g.Z = tX * tY;

#ifdef USE_MIXED_ADDITION
    g.to_special();
#endif

    return in;
}

std::ostream& operator<<(std::ostream& out, const std::vector<jubjub_bn254_G1> &v)
{
    out << v.size() << "\n";
    for (const jubjub_bn254_G1& t : v)
    {
        out << t << OUTPUT_NEWLINE;
    }

    return out;
}

std::istream& operator>>(std::istream& in, std::vector<jubjub_bn254_G1> &v)
{
    v.clear();

    size_t s;
    in >> s;
    v.reserve(s);
    consume_newline(in);

    for (size_t i = 0; i < s; ++i)
    {
        jubjub_bn254_G1 g;
        in >> g;
        v.emplace_back(g);
        consume_OUTPUT_NEWLINE(in);
    }

    return in;
}

void jubjub_bn254_G1::batch_to_special_all_non_zeros(std::vector<jubjub_bn254_G1> &vec)
{
    std::vector<jubjub_bn254_Fq> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el: vec)
    {
        Z_vec.emplace_back(el.Z);
    }
    batch_invert<jubjub_bn254_Fq>(Z_vec);

    const jubjub_bn254_Fq one = jubjub_bn254_Fq::one();

    for (size_t i = 0; i < vec.size(); ++i)
    {
        vec[i].X = vec[i].X * Z_vec[i];
        vec[i].Y = vec[i].Y * Z_vec[i];
        vec[i].Z = one;
    }
}

} // namespace libff