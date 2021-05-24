/*
 * MIT License
 *
 * Parametric Cubic Spline Library
 * Copyright (c) 2021-, Michael Heidingsfeld
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include <initializer_list>
#include <vector>
#include <cmath>

#include "gtest/gtest.h"
#include "parametric_cubic_spline/parametric_cubic_spline.h"

using namespace parametric_cubic_spline;

#define PIVOT_POINTS    {1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0,-1.0}
#define EVAL_POINTS     {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
#define EXPECTED_VALUES {1.0000, 0.0000, 0.1634,-0.1274,-0.5328,-0.1792, \
                        -0.9482,-0.0798,-0.9600, 0.2320,-0.6500, 0.6500, \
                        -0.2320, 0.9600, 0.0798, 0.9482, 0.1792, 0.5328, \
                         0.1274,-0.1634, 0.0000,-1.0000};

class TestProblem
{
public:
    std::vector<float> points_;
    std::size_t num_points_;
    std::size_t num_dims_;
    BoundaryCondition left_bc_;
    BoundaryCondition right_bc_;
    std::vector<float> left_tangent_;
    std::vector<float> right_tangent_;
    std::vector<float> eval_pos_;
    std::vector<float> expected_points_;

    TestProblem() = default;
    TestProblem(
        std::initializer_list<float> points, 
        std::size_t num_points,
        std::size_t num_dims,
        BoundaryCondition left_bc,
        BoundaryCondition right_bc,
        std::initializer_list<float> left_tangent,
        std::initializer_list<float> right_tangent,
        std::initializer_list<float> eval_pos,
        std::initializer_list<float> expected_points
    ) :
        points_(points),
        num_points_(num_points),
        num_dims_(num_dims),
        left_bc_(left_bc),            
        right_bc_(right_bc),
        left_tangent_(left_tangent),
        right_tangent_(right_tangent),
        eval_pos_(eval_pos),
        expected_points_(expected_points)
    {};   
};

static const TestProblem test_problem1(
    { 1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0,-1.0}, 
    4, 
    2,
    BoundaryCondition::Natural,
    BoundaryCondition::Natural,
    { 0.0, 0.0},
    { 0.0, 0.0},    
    { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 },
    { 1.0000, 0.0000, 0.1634,-0.1274,-0.5328,-0.1792, \
     -0.9482,-0.0798,-0.9600, 0.2320,-0.6500, 0.6500, \
     -0.2320, 0.9600, 0.0798, 0.9482, 0.1792, 0.5328, \
      0.1274,-0.1634, 0.0000,-1.0000 }
);
static const TestProblem test_problem2(
    { 1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0,-1.0}, 
    4, 
    2,
    BoundaryCondition::Hermite,
    BoundaryCondition::Hermite,
    { 0.0,-1.0},
    {-1.0, 0.0},
    { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 },
    { 1.0000, 0.0000, 0.6352,-0.2268,-0.1424,-0.2784, \
     -0.8576,-0.1116,-1.0731, 0.3003,-0.7917, 0.7917, \
     -0.3003, 1.0731, 0.1116, 0.8576, 0.2784, 0.1424, \
      0.2268,-0.6352, 0.0000,-1.0000 }
);


class MyTestFixture: public ::testing::TestWithParam<TestProblem> { 
public: 
   void SetUp() {}
   void TestBody() {}
   void TearDown() {}
};

TEST_P(MyTestFixture, MyTestName)
{
    TestProblem problem = GetParam();

    Spline<float, Dynamic, Dynamic> spline;
    spline.set(
        problem.points_.data(), 
        problem.num_points_, 
        problem.num_dims_,
        problem.left_bc_,
        problem.right_bc_,
        problem.left_tangent_.data(),
        problem.right_tangent_.data()
    );

    std::size_t eval_points_size = problem.eval_pos_.size()*problem.num_dims_;
    std::vector<float> eval_points(eval_points_size, 0.0);
    spline.eval(problem.eval_pos_.data(), 11, eval_points.data());        

    for(std::size_t i = 0; i < eval_points_size; i++)
    {
        std::cout << eval_points[i] << " | " << problem.expected_points_[i] << std::endl;
        EXPECT_LT(fabs(eval_points[i] - problem.expected_points_[i]), 0.001);
    }
}

INSTANTIATE_TEST_SUITE_P(
    MyTestSuite,
    MyTestFixture,
    ::testing::Values(
        test_problem1,
        test_problem2
    )
);

/*
TEST(ParametricCubicSplineTests, VectorDynamicPointsDynamicDims)
{
    std::vector<double> p = PIVOT_POINTS;
    std::vector<double> s = EVAL_POINTS;
    std::vector<double> x_expected = EXPECTED_VALUES;
    std::vector<double> x(2*11, 0.0);

    Spline<double, Dynamic, Dynamic> spline;
    spline.set(p.data(), 4, 2);
    spline.eval(s.data(), 11, x.data());

    for(std::size_t i = 0; i < x.size(); i++)
    {
        EXPECT_LT(fabs(x[i] - x_expected[i]), 0.001);
    }
}

TEST(ParametricCubicSplineTests, VectorDynamicPointsDynamicDimsHermite)
{
    std::vector<double> p = PIVOT_POINTS;
    std::vector<double> s = EVAL_POINTS;
    std::vector<double> x_expected = EXPECTED_VALUES;
    std::vector<double> x(2*11, 0.0);
    std::vector<double> zero_bc = {0.0, 0.0};

    Spline<double, Dynamic, Dynamic> spline;
    spline.set(p.data(), 4, 2, BoundaryCondition::Hermite, BoundaryCondition::Hermite, zero_bc.data(), zero_bc.data());
    spline.eval(s.data(), 11, x.data());

    for(std::size_t i = 0; i < x.size(); i++)
    {
        EXPECT_LT(fabs(x[i] - x_expected[i]), 0.001);
    }
}

TEST(ParametricCubicSplineTests, VectorDynamicPointsFixedDims)
{
    std::vector<double> p = PIVOT_POINTS;
    std::vector<double> s = EVAL_POINTS;
    std::vector<double> x_expected = EXPECTED_VALUES;
    std::vector<double> x(2*11, 0.0);

    Spline<double, Dynamic, 2> spline;
    spline.set(p.data(), 4);
    spline.eval(s.data(), 11, x.data());

    for(std::size_t i = 0; i < 2*11; i++)
    {
        EXPECT_LT(fabs(x[i] - x_expected[i]), 0.001);
    }
}

TEST(ParametricCubicSplineTests, DynamicCArray)
{
    double p[2*4] = PIVOT_POINTS;
    double s[11] = EVAL_POINTS;
    double x_expected[2*11] = EXPECTED_VALUES;
    double x[2*11] = {0.0};

    Spline<double, Dynamic, Dynamic> spline;
    spline.set(p, 4, 2);
    spline.eval(s, 11, x);

    for(std::size_t i = 0; i < 2*11; i++)
    {
        EXPECT_LT(fabs(x[i] - x_expected[i]), 0.001);
    }
}

TEST(ParametricCubicSplineTests, StaticCArray)
{
    double p[2*4] = PIVOT_POINTS;
    double s[11] = EVAL_POINTS;
    double x_expected[2*11] = EXPECTED_VALUES;
    double x[2*11] = {0.0};

    Spline<double, 4, 2> spline;
    spline.set(p);
    spline.eval(s, 11, x);

    for(std::size_t i = 0; i < 2*11; i++)
    {
        EXPECT_LT(fabs(x[i] - x_expected[i]), 0.001);
    }
}
*/