#ifndef JUBJUB_ED181_PP_HPP_
#define JUBJUB_ED181_PP_HPP_
#include <libff/algebra/curves/jubjub_ed181/jubjub_ed181_g1.hpp>
#include <libff/algebra/curves/jubjub_ed181/jubjub_ed181_init.hpp>
#include <libff/algebra/curves/edwards181/edwards181_fields.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class jubjub_ed181_pp {
public:
    typedef jubjub_ed181_Fr Fp_type;
    typedef jubjub_ed181_G1 G1_type;
    typedef jubjub_ed181_Fq Fq_type;

    static void init_public_params();
    static edwards181_Fr inner2outer(const jubjub_ed181_Fq &x);
    static jubjub_ed181_Fq outer2inner(const edwards181_Fr &x);
};

} // libff

#endif // JUBJUB_ED181_PP_HPP_
