#include <libff/algebra/curves/jubjub_ed97/jubjub_ed97_pp.hpp>

namespace libff {

void jubjub_ed97_pp::init_public_params()
{
    init_jubjub_ed97_params();
}

edwards97_Fr jubjub_ed97_pp::inner2outer(const jubjub_ed97_Fq &x){
    edwards97_Fr xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

jubjub_ed97_Fq jubjub_ed97_pp::outer2inner(const edwards97_Fr &x){
    jubjub_ed97_Fq xo;

    xo.mont_repr = x.mont_repr;
    return xo;
}

} // libff
