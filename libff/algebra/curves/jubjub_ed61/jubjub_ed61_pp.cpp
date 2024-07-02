#include <libff/algebra/curves/jubjub_ed61/jubjub_ed61_pp.hpp>

namespace libff {

void jubjub_ed61_pp::init_public_params()
{
    init_jubjub_ed61_params();
}

edwards61_Fr jubjub_ed61_pp::inner2outer(const jubjub_ed61_Fq &x){
    edwards61_Fr xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

jubjub_ed61_Fq jubjub_ed61_pp::outer2inner(const edwards61_Fr &x){
    jubjub_ed61_Fq xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

} // libff
