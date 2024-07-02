#include <libff/algebra/curves/baby_jubjub/baby_jubjub_pp.hpp>

namespace libff {

void baby_jubjub_pp::init_public_params()
{
    init_baby_jubjub_params();
}

alt_bn128_Fr baby_jubjub_pp::inner2outer(const baby_jubjub_Fq &x){
    alt_bn128_Fr xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

baby_jubjub_Fq baby_jubjub_pp::outer2inner(const alt_bn128_Fr &x){
    baby_jubjub_Fq xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

} // libff
