#ifndef JUBJUB_ED97_PP_HPP_
#define JUBJUB_ED97_PP_HPP_
#include <libff/algebra/curves/jubjub_ed97/jubjub_ed97_g1.hpp>
#include <libff/algebra/curves/jubjub_ed97/jubjub_ed97_init.hpp>
#include <libff/algebra/curves/edwards97/edwards97_fields.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class jubjub_ed97_pp {
public:
    typedef jubjub_ed97_Fr Fp_type;
    typedef jubjub_ed97_G1 G1_type;
    typedef jubjub_ed97_Fq Fq_type;

    static void init_public_params();
    static edwards97_Fr inner2outer(const jubjub_ed97_Fq &x);
    static jubjub_ed97_Fq outer2inner(const edwards97_Fr &x);
};

} // libff

#endif // JUBJUB_ED97_PP_HPP_
