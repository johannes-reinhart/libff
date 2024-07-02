#include <libff/algebra/curves/jubjub_bn124/jubjub_bn124_pp.hpp>

namespace libff {

void jubjub_bn124_pp::init_public_params()
{
    init_jubjub_bn124_params();
}

bn124_Fr jubjub_bn124_pp::inner2outer(const jubjub_bn124_Fq &x){
    bn124_Fr xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

jubjub_bn124_Fq jubjub_bn124_pp::outer2inner(const bn124_Fr &x){
    jubjub_bn124_Fq xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

} // libff
