/**
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <iostream>
#include <gtest/gtest.h>

#ifdef MULTICORE
#include <omp.h>
#endif

#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>

#include <libff/algebra/scalar_multiplication/multiexp.hpp>

using namespace libff;

using std::size_t;

class MultiexpTest: public ::testing::Test {
public:
    MultiexpTest()
    {
        alt_bn128_pp::init_public_params();
    }
};


template<typename GroupT, typename FieldT>
void test_multiexp_plain()
{
#ifdef MULTICORE
    const size_t chunks = omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call omp_set_num_threads()
#else
    const size_t chunks = 1;
#endif

    GroupT result_naive = GroupT::zero();
    GroupT result_naive_plain = GroupT::zero();
    GroupT result_bos_coster = GroupT::zero();
    GroupT result_BDLO12 = GroupT::zero();
    GroupT result_expected = GroupT::zero();

    FieldT r1 = FieldT("76749407");
    FieldT r2 = FieldT("44410867");
    FieldT r3 = FieldT("53490921");

    FieldT s1 = FieldT("1234");
    FieldT s2 = FieldT("4056764");
    FieldT s3 = FieldT("1945");

    std::vector<GroupT> values;
    std::vector<FieldT> scalars;

    values.push_back(r1 * GroupT::one());
    values.push_back(r2 * GroupT::one());
    values.push_back(r3 * GroupT::one());

    scalars.push_back(s1);
    scalars.push_back(s2);
    scalars.push_back(s3);

    result_naive = multi_exp<GroupT, FieldT, libff::multi_exp_method_naive_plain>(
                values.begin(),
                values.end(),
                scalars.begin(),
                scalars.end(),
                chunks);

    result_naive_plain = multi_exp<GroupT, FieldT, libff::multi_exp_method_naive_plain>(
                    values.begin(),
                    values.end(),
                    scalars.begin(),
                    scalars.end(),
                    chunks);

    result_bos_coster = multi_exp<GroupT, FieldT, libff::multi_exp_method_bos_coster>(
                    values.begin(),
                    values.end(),
                    scalars.begin(),
                    scalars.end(),
                    chunks);

    result_BDLO12 = multi_exp<GroupT, FieldT, libff::multi_exp_method_BDLO12>(
                values.begin(),
                values.end(),
                scalars.begin(),
                scalars.end(),
                chunks);

    result_expected = (s1*r1 + s2*r2 + s3*r3) * GroupT::one();

    EXPECT_EQ(result_naive, result_expected);
    EXPECT_EQ(result_naive_plain, result_expected);
    EXPECT_EQ(result_bos_coster, result_expected);
    EXPECT_EQ(result_BDLO12, result_expected);
}

TEST_F(MultiexpTest, PlainTest)
{
    test_multiexp_plain<G1<alt_bn128_pp>, Fr<alt_bn128_pp>>();
}
