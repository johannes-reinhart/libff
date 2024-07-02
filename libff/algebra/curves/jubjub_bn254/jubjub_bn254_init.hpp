#ifndef JUBJUB_BN254_INIT_HPP_
#define JUBJUB_BN254_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/jubjub_bn254/jubjub_bn254_fields.hpp>

namespace libff {

void init_jubjub_bn254_params();

class jubjub_bn254_G1;

} // libff
#endif // JUBJUB_BN254_INIT_HPP_
