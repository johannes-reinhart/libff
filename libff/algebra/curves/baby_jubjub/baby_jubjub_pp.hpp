#ifndef BABY_JUBJUB_PP_HPP_
#define BABY_JUBJUB_PP_HPP_
#include <libff/algebra/curves/baby_jubjub/baby_jubjub_g1.hpp>
#include <libff/algebra/curves/baby_jubjub/baby_jubjub_init.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_fields.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class baby_jubjub_pp {
public:
    typedef baby_jubjub_Fr Fp_type;
    typedef baby_jubjub_G1 G1_type;
    typedef baby_jubjub_Fq Fq_type;

    static void init_public_params();
    static alt_bn128_Fr inner2outer(const baby_jubjub_Fq &x);
    static baby_jubjub_Fq outer2inner(const alt_bn128_Fr &x);
};

} // libff

#endif // BABY_JUBJUB_PP_HPP_
