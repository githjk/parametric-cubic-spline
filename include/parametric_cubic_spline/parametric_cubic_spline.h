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
#pragma once

namespace parametric_cubic_spline {

namespace internal {

    template<typename T, std::size_t N>
    class StorageType;

} // namespace: internal

/**
 * Constant used to express dynamic size
 */
static const std::size_t Dynamic = 0;

/**
 * Boundary condition class
 */
enum class BoundaryCondition
{
    Natural,
    Hermite,
    Periodic,
    NotAKnot
};

/**
 * Spline class
 */
template<
    typename T,
    std::size_t NumPoints = Dynamic,
    std::size_t NumDims = Dynamic
>
class Spline
{
    std::size_t num_points_;
    std::size_t num_dims_;
    const T *points_;
    internal::StorageType<T, NumPoints*NumDims> moments_;

public:
    Spline();

    // variable points, variable dims, optional bc
    void set(
        const T *points,
        const std::size_t num_points,
        const std::size_t num_dims,        
        const BoundaryCondition left_bc = BoundaryCondition::Natural,
        const BoundaryCondition right_bc = BoundaryCondition::Natural,
        const T *left_tangent = nullptr,
        const T *right_tangent = nullptr
    );

    // variable points, fixed dims, optional bc
    void set(
        const T *points,
        const std::size_t num_points,
        const BoundaryCondition left_bc = BoundaryCondition::Natural,
        const BoundaryCondition right_bc = BoundaryCondition::Natural,
        const T *left_tangent = nullptr,
        const T *right_tangent = nullptr
    );

    // fixed points, fixed dims, optional bc
    void set(
        const T *points,
        const BoundaryCondition left_bc = BoundaryCondition::Natural,
        const BoundaryCondition right_bc = BoundaryCondition::Natural,
        const T *left_tangent = nullptr,
        const T *right_tangent = nullptr
    );

    // variable lengths
    void eval(
        const T *pos,
        const std::size_t num_pos,        
        T *out_points
    );

    // single point
    void eval(
        const T pos,
        T *out_point
    );

private:
    static void compute_moments(
        const T* points,
        const std::size_t num_points,
        const std::size_t num_dims,        
        const BoundaryCondition left_bc,
        const BoundaryCondition right_bc,
        const T* left_tangent,
        const T* right_tangent,
        internal::StorageType<T, NumPoints*NumDims> &m
    );

    static void tdma(
        const std::size_t num_points,
        const std::size_t num_dims,
        internal::StorageType<T, NumPoints> &a,
        internal::StorageType<T, NumPoints> &b,
        internal::StorageType<T, NumPoints> &c,
        internal::StorageType<T, NumPoints*NumDims> &d,
        T *q = nullptr,
        T *u = nullptr      
    );
};

} // namespace: parametric_cubic_spline

#include "parametric_cubic_spline/impl/parametric_cublic_spline.hpp"