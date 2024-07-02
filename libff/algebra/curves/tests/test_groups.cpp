/**
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <iostream>
#include <typeinfo>

#include <gtest/gtest.h>

#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/algebra/curves/edwards58/edwards58_pp.hpp>
#include <libff/algebra/curves/edwards61/edwards61_pp.hpp>
#include <libff/algebra/curves/edwards97/edwards97_pp.hpp>
#include <libff/algebra/curves/edwards181/edwards181_pp.hpp>

#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/curves/bn124/bn124_pp.hpp>
#include <libff/algebra/curves/bn183/bn183_pp.hpp>
#include <libff/algebra/curves/bn254/bn254_pp.hpp>

#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_pp.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>

#include <libff/algebra/curves/baby_jubjub/baby_jubjub_pp.hpp>
#include <libff/algebra/curves/jubjub_bn124/jubjub_bn124_pp.hpp>
#include <libff/algebra/curves/jubjub_bn183/jubjub_bn183_pp.hpp>
#include <libff/algebra/curves/jubjub_bn254/jubjub_bn254_pp.hpp>
#include <libff/algebra/curves/jubjub_ed58/jubjub_ed58_pp.hpp>
#include <libff/algebra/curves/jubjub_ed61/jubjub_ed61_pp.hpp>
#include <libff/algebra/curves/jubjub_ed97/jubjub_ed97_pp.hpp>
#include <libff/algebra/curves/jubjub_ed181/jubjub_ed181_pp.hpp>


#ifdef CURVE_BN128
#include <libff/algebra/curves/bn128/bn128_pp.hpp>
#endif
#include <sstream>


using namespace libff;

using std::size_t;

class CurveGroupsTest: public ::testing::Test {
public:
    CurveGroupsTest()
    {
        edwards_pp::init_public_params();
        edwards58_pp::init_public_params();
        edwards61_pp::init_public_params();
        edwards97_pp::init_public_params();
        edwards181_pp::init_public_params();

        mnt4_pp::init_public_params();
        mnt6_pp::init_public_params();

        alt_bn128_pp::init_public_params();
        bn124_pp::init_public_params();
        bn183_pp::init_public_params();
        bn254_pp::init_public_params();

        bls12_381_pp::init_public_params();

        baby_jubjub_pp::init_public_params();
        jubjub_bn124_pp::init_public_params();
        jubjub_bn183_pp::init_public_params();
        jubjub_bn254_pp::init_public_params();
        jubjub_ed58_pp::init_public_params();
        jubjub_ed61_pp::init_public_params();
        jubjub_ed97_pp::init_public_params();
        jubjub_ed181_pp::init_public_params();

#ifdef CURVE_BN128 // BN128 has fancy dependencies so it may be disabled
        bn128_pp::init_public_params();
#endif
    }
};

template<typename GroupT>
void test_mixed_add()
{
    GroupT base, el, result;

    base = GroupT::zero();
    el = GroupT::zero();
    el.to_special();
    result = base.mixed_add(el);
    EXPECT_EQ(result, base + el);

    base = GroupT::zero();
    el = GroupT::random_element();
    el.to_special();
    result = base.mixed_add(el);
    EXPECT_EQ(result, base + el);

    base = GroupT::random_element();
    el = GroupT::zero();
    el.to_special();
    result = base.mixed_add(el);
    EXPECT_EQ(result, base + el);

    base = GroupT::random_element();
    el = GroupT::random_element();
    el.to_special();
    result = base.mixed_add(el);
    EXPECT_EQ(result, base + el);

    base = GroupT::random_element();
    el = base;
    el.to_special();
    result = base.mixed_add(el);
    EXPECT_EQ(result, base.dbl());
}

template<typename GroupT>
void test_group()
{
    bigint<1> rand1 = bigint<1>("76749407");
    bigint<1> rand2 = bigint<1>("44410867");
    bigint<1> randsum = bigint<1>("121160274");

    GroupT zero = GroupT::zero();
    EXPECT_EQ(zero, zero);
    GroupT one = GroupT::one();
    EXPECT_EQ(one, one);
    GroupT two = bigint<1>(2L) * GroupT::one();
    EXPECT_EQ(two, two);
    GroupT five = bigint<1>(5L) * GroupT::one();

    GroupT three = bigint<1>(3L) * GroupT::one();
    GroupT four = bigint<1>(4L) * GroupT::one();

    EXPECT_EQ(two+five, three+four);

    GroupT a = random_element_non_zero_one<GroupT>();
    GroupT b = random_element_non_zero_one<GroupT>();

    EXPECT_NE(one, zero);
    EXPECT_NE(a, zero);
    EXPECT_NE(a, one);
    EXPECT_NE(b, zero);
    EXPECT_NE(b, one);

    EXPECT_EQ(a.dbl(), a + a);
    EXPECT_EQ(b.dbl(), b + b);
    EXPECT_EQ(one.add(two), three);
    EXPECT_EQ(two.add(one), three);
    EXPECT_EQ(a + b, b + a);
    EXPECT_EQ(a - a, zero);
    EXPECT_EQ(a - b, a + (-b));
    EXPECT_EQ(a - b, (-b) + a);

    // handle special cases
    EXPECT_EQ(zero + (-a), -a);
    EXPECT_EQ(zero - a, -a);
    EXPECT_EQ(a - zero, a);
    EXPECT_EQ(a + zero, a);
    EXPECT_EQ(zero + a, a);

    EXPECT_EQ((a + b).dbl(), (a + b) + (b + a));
    EXPECT_EQ(bigint<1>("2") * (a + b), (a + b) + (b + a));

    EXPECT_EQ(rand1 * a + rand2 * a, randsum * a);

    EXPECT_EQ(GroupT::order() * a, zero);
    EXPECT_EQ(GroupT::order() * one, zero);
    EXPECT_NE(GroupT::order() * a - a, zero);
    EXPECT_NE(GroupT::order() * one - one, zero);

    test_mixed_add<GroupT>();
}

template<typename GroupT>
void test_mul_by_q()
{
    GroupT a = GroupT::random_element();
    EXPECT_EQ(GroupT::field_char() * a, a.mul_by_q()) << "Type of a: " << typeid(a).name();
}

template<typename GroupT>
void test_output()
{
    GroupT g = GroupT::zero();

    /* ate-pairing contained optimizations specific to the original curve that were breaking
       point addition with extremely small probability, so this code was run for 1000 times
       in case there was a missing carry. Since no problems were found, this is now reduced
       to only 10 times for quick testing. */
    for (size_t i = 0; i < 10; ++i)
    {
        GroupT g_ser = reserialize(g);
        EXPECT_EQ(g, g_ser);
        // Use a random point in next iteration
        g = GroupT::random_element();
    }
}

TEST_F(CurveGroupsTest, GroupTest)
{
    test_group<G1<edwards_pp> >();
    test_group<G2<edwards_pp> >();

    test_group<G1<edwards58_pp> >();
    test_group<G2<edwards58_pp> >();

    test_group<G1<edwards61_pp> >();
    test_group<G2<edwards61_pp> >();

    test_group<G1<edwards97_pp> >();
    test_group<G2<edwards97_pp> >();

    test_group<G1<edwards181_pp> >();
    test_group<G2<edwards181_pp> >();

    test_group<G1<mnt4_pp> >();
    test_group<G2<mnt4_pp> >();

    test_group<G1<mnt6_pp> >();
    test_group<G2<mnt6_pp> >();

    test_group<G1<alt_bn128_pp> >();
    test_group<G2<alt_bn128_pp> >();
    test_group<G1<bn124_pp> >();
    test_group<G2<bn124_pp> >();
    test_group<G1<bn183_pp> >();
    test_group<G2<bn183_pp> >();
    test_group<G1<bn254_pp> >();
    test_group<G2<bn254_pp> >();

    test_group<G1<bls12_381_pp> >();
    test_group<G2<bls12_381_pp> >();

    test_group<G1<baby_jubjub_pp> >();
    test_group<G1<jubjub_bn124_pp> >();
    test_group<G1<jubjub_bn183_pp> >();
    test_group<G1<jubjub_bn254_pp> >();
    test_group<G1<jubjub_ed58_pp> >();
    test_group<G1<jubjub_ed61_pp> >();
    test_group<G1<jubjub_ed97_pp> >();
    test_group<G1<jubjub_ed181_pp> >();

#ifdef CURVE_BN128       // BN128 has fancy dependencies so it may be disabled
    test_group<G1<bn128_pp> >();
    test_group<G2<bn128_pp> >();
#endif
}

TEST_F(CurveGroupsTest, OutputTest)
{
    test_output<G1<edwards_pp> >();
    test_output<G2<edwards_pp> >();

    test_output<G1<edwards58_pp> >();
    test_output<G2<edwards58_pp> >();

    test_output<G1<edwards61_pp> >();
    test_output<G2<edwards61_pp> >();

    test_output<G1<edwards97_pp> >();
    test_output<G2<edwards97_pp> >();

    test_output<G1<edwards181_pp> >();
    test_output<G2<edwards181_pp> >();

    test_output<G1<mnt4_pp> >();
    test_output<G2<mnt4_pp> >();

    test_output<G1<mnt6_pp> >();
    test_output<G2<mnt6_pp> >();

    test_output<G1<alt_bn128_pp> >();
    test_output<G2<alt_bn128_pp> >();
    test_output<G1<bn124_pp> >();
    test_output<G2<bn124_pp> >();
    test_output<G1<bn183_pp> >();
    test_output<G2<bn183_pp> >();
    test_output<G1<bn254_pp> >();
    test_output<G2<bn254_pp> >();

    test_output<G1<bls12_381_pp> >();
    test_output<G2<bls12_381_pp> >();

    test_output<G1<baby_jubjub_pp> >();
    test_output<G1<jubjub_bn124_pp> >();
    test_output<G1<jubjub_bn183_pp> >();
    test_output<G1<jubjub_bn254_pp> >();
    test_output<G1<jubjub_ed58_pp> >();
    test_output<G1<jubjub_ed61_pp> >();
    test_output<G1<jubjub_ed97_pp> >();
    test_output<G1<jubjub_ed181_pp> >();

#ifdef CURVE_BN128       // BN128 has fancy dependencies so it may be disabled
    test_output<G1<bn128_pp> >();
    test_output<G2<bn128_pp> >();
#endif
}

TEST_F(CurveGroupsTest, MulByQTest)
{
    test_mul_by_q<G2<edwards_pp> >();
    test_mul_by_q<G2<edwards58_pp> >();
    test_mul_by_q<G2<edwards61_pp> >();
    test_mul_by_q<G2<edwards97_pp> >();
    test_mul_by_q<G2<edwards181_pp> >();

    test_mul_by_q<G2<mnt4_pp> >();
    test_mul_by_q<G2<mnt6_pp> >();
    test_mul_by_q<G2<alt_bn128_pp> >();
    test_mul_by_q<G2<bn124_pp> >();
    test_mul_by_q<G2<bn183_pp> >();
    test_mul_by_q<G2<bn254_pp> >();

    test_mul_by_q<G2<bls12_381_pp> >();
}
