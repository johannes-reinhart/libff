#ifndef JUBJUB_ED58_PP_HPP_
#define JUBJUB_ED58_PP_HPP_
#include <libff/algebra/curves/jubjub_ed58/jubjub_ed58_g1.hpp>
#include <libff/algebra/curves/jubjub_ed58/jubjub_ed58_init.hpp>
#include <libff/algebra/curves/edwards58/edwards58_fields.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class jubjub_ed58_pp {
public:
    typedef jubjub_ed58_Fr Fp_type;
    typedef jubjub_ed58_G1 G1_type;
    typedef jubjub_ed58_Fq Fq_type;

    static void init_public_params();
    static edwards58_Fr inner2outer(const jubjub_ed58_Fq &x);
    static jubjub_ed58_Fq outer2inner(const edwards58_Fr &x);
};

} // libff

#endif // JUBJUB_ED58_PP_HPP_
