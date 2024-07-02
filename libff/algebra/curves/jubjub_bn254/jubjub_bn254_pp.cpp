#include <libff/algebra/curves/jubjub_bn254/jubjub_bn254_pp.hpp>

namespace libff {

void jubjub_bn254_pp::init_public_params()
{
    init_jubjub_bn254_params();
}

bn254_Fr jubjub_bn254_pp::inner2outer(const jubjub_bn254_Fq &x){
    bn254_Fr xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

jubjub_bn254_Fq jubjub_bn254_pp::outer2inner(const bn254_Fr &x){
    jubjub_bn254_Fq xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

} // libff
