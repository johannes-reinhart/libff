#include <libff/algebra/curves/jubjub_bn183/jubjub_bn183_pp.hpp>

namespace libff {

void jubjub_bn183_pp::init_public_params()
{
    init_jubjub_bn183_params();
}

bn183_Fr jubjub_bn183_pp::inner2outer(const jubjub_bn183_Fq &x){
    bn183_Fr xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

jubjub_bn183_Fq jubjub_bn183_pp::outer2inner(const bn183_Fr &x){
    jubjub_bn183_Fq xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

} // libff
