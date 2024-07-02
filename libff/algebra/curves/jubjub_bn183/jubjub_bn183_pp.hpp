#ifndef JUBJUB_BN183_PP_HPP_
#define JUBJUB_BN183_PP_HPP_
#include <libff/algebra/curves/jubjub_bn183/jubjub_bn183_g1.hpp>
#include <libff/algebra/curves/jubjub_bn183/jubjub_bn183_init.hpp>
#include <libff/algebra/curves/bn183/bn183_fields.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class jubjub_bn183_pp {
public:
    typedef jubjub_bn183_Fr Fp_type;
    typedef jubjub_bn183_G1 G1_type;
    typedef jubjub_bn183_Fq Fq_type;

    static void init_public_params();
    static bn183_Fr inner2outer(const jubjub_bn183_Fq &x);
    static jubjub_bn183_Fq outer2inner(const bn183_Fr &x);
};

} // libff

#endif // JUBJUB_BN183_PP_HPP_
