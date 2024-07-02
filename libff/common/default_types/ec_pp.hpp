/** @file
 *****************************************************************************
 This file defines default_ec_pp based on the CURVE=... make flag, which selects
 which elliptic curve is used to implement group arithmetic and pairings.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef EC_PP_HPP_
#define EC_PP_HPP_

/************************ Pick the elliptic curve ****************************/

#include <libff/common/default_types/ec_aliases.hpp>

#ifdef CURVE_BLS12_381
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
namespace libff {
typedef bls12_381_pp default_ec_pp;
} // namespace libff
#endif

#ifdef CURVE_ALT_BN128
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
namespace libff {
typedef alt_bn128_pp default_ec_pp;
} // namespace libff
#endif

#ifdef CURVE_BN128
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/bn128/bn128_pp.hpp>
namespace libff {
typedef bn128_pp default_ec_pp;
} // namespace libff
#endif

#ifdef CURVE_BN124
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/bn124/bn124_pp.hpp>
namespace libff {
typedef bn124_pp default_ec_pp;
} // namespace libff
#endif

#ifdef CURVE_BN183
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/bn183/bn183_pp.hpp>
namespace libff {
typedef bn183_pp default_ec_pp;
} // namespace libff
#endif

#ifdef CURVE_BN254
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/bn254/bn254_pp.hpp>
namespace libff {
typedef bn254_pp default_ec_pp;
} // namespace libff
#endif

#ifdef CURVE_EDWARDS
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/edwards/edwards_pp.hpp>
namespace libff {
typedef edwards_pp default_ec_pp;
} // namespace libff
#endif

#ifdef CURVE_EDWARDS58
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/edwards58/edwards58_pp.hpp>
namespace libff {
    typedef edwards58_pp default_ec_pp;
} // namespace libff
#endif

#ifdef CURVE_EDWARDS61
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/edwards61/edwards61_pp.hpp>
namespace libff {
    typedef edwards61_pp default_ec_pp;
} // namespace libff
#endif

#ifdef CURVE_EDWARDS97
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/edwards97/edwards97_pp.hpp>
namespace libff {
    typedef edwards97_pp default_ec_pp;
} // namespace libff
#endif

#ifdef CURVE_EDWARDS181
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/edwards181/edwards181_pp.hpp>
namespace libff {
typedef edwards181_pp default_ec_pp;
} // namespace libff
#endif



#ifdef CURVE_MNT4
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/mnt/mnt4/mnt4_pp.hpp>
namespace libff {
typedef mnt4_pp default_ec_pp;
} // namespace libff
#endif

#ifdef CURVE_MNT6
#define LIBFF_DEFAULT_EC_PP_DEFINED
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>
namespace libff {
typedef mnt6_pp default_ec_pp;
} // namespace libff
#endif

#ifndef LIBFF_DEFAULT_EC_PP_DEFINED
#error You must define one of the CURVE_* symbols to pick a curve for pairings.
#endif

#endif // EC_PP_HPP_
