#ifndef JUBJUB_BN254_PP_HPP_
#define JUBJUB_BN254_PP_HPP_
#include <libff/algebra/curves/jubjub_bn254/jubjub_bn254_g1.hpp>
#include <libff/algebra/curves/jubjub_bn254/jubjub_bn254_init.hpp>
#include <libff/algebra/curves/bn254/bn254_fields.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class jubjub_bn254_pp {
public:
    typedef jubjub_bn254_Fr Fp_type;
    typedef jubjub_bn254_G1 G1_type;
    typedef jubjub_bn254_Fq Fq_type;

    static void init_public_params();
    static bn254_Fr inner2outer(const jubjub_bn254_Fq &x);
    static jubjub_bn254_Fq outer2inner(const bn254_Fr &x);
};

} // libff

#endif // JUBJUB_BN254_PP_HPP_
