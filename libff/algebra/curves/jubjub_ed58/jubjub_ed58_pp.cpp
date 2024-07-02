#include <libff/algebra/curves/jubjub_ed58/jubjub_ed58_pp.hpp>

namespace libff {

void jubjub_ed58_pp::init_public_params()
{
    init_jubjub_ed58_params();
}

edwards58_Fr jubjub_ed58_pp::inner2outer(const jubjub_ed58_Fq &x){
    edwards58_Fr xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

jubjub_ed58_Fq jubjub_ed58_pp::outer2inner(const edwards58_Fr &x){
    jubjub_ed58_Fq xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

} // libff
