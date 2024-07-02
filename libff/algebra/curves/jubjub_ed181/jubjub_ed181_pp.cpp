#include <libff/algebra/curves/jubjub_ed181/jubjub_ed181_pp.hpp>

namespace libff {

void jubjub_ed181_pp::init_public_params()
{
    init_jubjub_ed181_params();
}

edwards181_Fr jubjub_ed181_pp::inner2outer(const jubjub_ed181_Fq &x){
    edwards181_Fr xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

jubjub_ed181_Fq jubjub_ed181_pp::outer2inner(const edwards181_Fr &x){
    jubjub_ed181_Fq xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

} // libff
