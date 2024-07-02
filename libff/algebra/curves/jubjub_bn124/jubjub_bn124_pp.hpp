#ifndef JUBJUB_BN124_PP_HPP_
#define JUBJUB_BN124_PP_HPP_
#include <libff/algebra/curves/jubjub_bn124/jubjub_bn124_g1.hpp>
#include <libff/algebra/curves/jubjub_bn124/jubjub_bn124_init.hpp>
#include <libff/algebra/curves/bn124/bn124_fields.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class jubjub_bn124_pp {
public:
    typedef jubjub_bn124_Fr Fp_type;
    typedef jubjub_bn124_G1 G1_type;
    typedef jubjub_bn124_Fq Fq_type;

    static void init_public_params();
    static bn124_Fr inner2outer(const jubjub_bn124_Fq &x);
    static jubjub_bn124_Fq outer2inner(const bn124_Fr &x);
};

} // libff

#endif // JUBJUB_BN124_PP_HPP_
