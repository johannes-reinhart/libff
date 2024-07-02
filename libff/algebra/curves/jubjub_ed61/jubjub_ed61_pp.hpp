#ifndef JUBJUB_ED61_PP_HPP_
#define JUBJUB_ED61_PP_HPP_
#include <libff/algebra/curves/jubjub_ed61/jubjub_ed61_g1.hpp>
#include <libff/algebra/curves/jubjub_ed61/jubjub_ed61_init.hpp>
#include <libff/algebra/curves/edwards61/edwards61_fields.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class jubjub_ed61_pp {
public:
    typedef jubjub_ed61_Fr Fp_type;
    typedef jubjub_ed61_G1 G1_type;
    typedef jubjub_ed61_Fq Fq_type;

    static void init_public_params();
    static edwards61_Fr inner2outer(const jubjub_ed61_Fq &x);
    static jubjub_ed61_Fq outer2inner(const edwards61_Fr &x);
};

} // libff

#endif // JUBJUB_ED61_PP_HPP_
