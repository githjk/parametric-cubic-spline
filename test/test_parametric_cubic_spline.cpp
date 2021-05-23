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