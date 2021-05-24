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

#include <cassert>

namespace parametric_cubic_spline {

namespace internal {

    /**
     * Statically allocated array
     */
    template<typename T, std::size_t N>
    class StorageType
    {
        std::array<T, N> data_;
    public:
        StorageType() = default;
        StorageType(std::size_t) { /* Do nothing */ }
        inline void resize(std::size_t) { /* Do nothing */ }
        inline T* data() { return data_.data(); }
        inline T& operator[](int pos) { return data_[pos]; }
        inline const T& operator[](int pos) const { return data_[pos]; }
    };

    /**
     * Partial template specialization for dynamically allocated array
     */
    template<typename T>
    class StorageType<T, Dynamic>
    {
        std::vector<T> data_;
    public:
        StorageType() = default;
        StorageType(std::size_t size) { resize(size); }
        inline void resize(std::size_t size) { data_ = std::vector<T>(size, 0.0); }
        inline T* data() { return data_.data(); }
        inline T& operator[](int pos) { return data_[pos]; }
        inline const T& operator[](int pos) const { return data_[pos]; }
    };

} // namespace: internal


template<typename T, std::size_t NumPoints, std::size_t NumDims>
Spline<T, NumPoints, NumDims>::Spline()
{
    static_assert(std::is_floating_point<T>::value, "T must be a floating point type.");
    static_assert(NumPoints != 1, "NumPoints must be either 'Dynamic' or greater than 1.");
}

template<typename T, std::size_t NumPoints, std::size_t NumDims>
void Spline<T, NumPoints, NumDims>::set(
    const T *points,
    const std::size_t num_points,
    const std::size_t num_dims,
    const BoundaryCondition left_bc,
    const BoundaryCondition right_bc,
    const T *left_tangent,
    const T *right_tangent
) {
    // Assign pointer to pivot points
    num_points_ = num_points;
    num_dims_ = num_dims;
    points_ = points;

    // In case of dynamic size, resize moments
    if(NumPoints == Dynamic || NumDims == Dynamic)
    {
        moments_.resize(num_points_*num_dims_);
    }

    // Compute moments
    compute_moments(points_, num_points_, num_dims_, left_bc, right_bc,
        left_tangent, right_tangent, moments_);
}

template<typename T, std::size_t NumPoints, std::size_t NumDims>
void Spline<T, NumPoints, NumDims>::set(
    const T *points,
    const std::size_t num_points,
    const BoundaryCondition left_bc,
    const BoundaryCondition right_bc,
    const T *left_tangent,
    const T *right_tangent
) {
    static_assert(NumDims > 0, "Number of dimensions 'NumDims' must be greater than zero.");

    set(points, num_points, NumDims, left_bc, right_bc, left_tangent, right_tangent);
}

template<typename T, std::size_t NumPoints, std::size_t NumDims>
void Spline<T, NumPoints, NumDims>::set(
    const T *points,
    const BoundaryCondition left_bc,
    const BoundaryCondition right_bc,
    const T *left_tangent,
    const T *right_tangent
) {
    static_assert(NumPoints > 1, "Number of points 'NumPoints' must be greater than one.");
    static_assert(NumDims > 0, "Number of dimensions 'NumDims' must be greater than zero.");

    set(points, NumPoints, NumDims, left_bc, right_bc, left_tangent, right_tangent);
}

template<typename T, std::size_t NumPoints, std::size_t NumDims>
void Spline<T, NumPoints, NumDims>::eval(
    const T pos,
    T *out_point
)
{
    std::size_t i = floor(pos * (num_points_ - 1));
    T t = fmod(pos * (num_points_ - 1), 1.0);
    if(i == num_points_ - 1)
    {
        i--;
        t = 1.0;
    }

    T t0 = pow(t, 3.0);
    T t1 = pow(1-t, 3.0);
    for(std::size_t j = 0; j < num_dims_; j++)
    {
        T c = (points_[(i+1)*num_dims_+j] - points_[i*num_dims_+j]) - 1.0/6.0*(moments_[(i+1)*num_dims_+j] - moments_[i*num_dims_+j]);
        T d = points_[i*num_dims_+j] - 1.0/6.0*moments_[i*num_dims_+j];
        *out_point++ = 1.0/6.0*(t1*moments_[i*num_dims_+j] + t0*moments_[(i+1)*num_dims_+j]) + c*t + d;
    }
}

template<typename T, std::size_t NumPoints, std::size_t NumDims>
void Spline<T, NumPoints, NumDims>::eval(
    const T *pos,
    const std::size_t num_pos,
    T *out_points
)
{
    for(std::size_t i = 0; i < num_pos; i++)
    {
        eval(pos[i], &(out_points[i*num_dims_]));
    }
}

template<typename T, std::size_t NumPoints, std::size_t NumDims>
void Spline<T, NumPoints, NumDims>::compute_moments(
    const T* points,
    const std::size_t num_points,
    const std::size_t num_dims,
    const BoundaryCondition left_bc,
    const BoundaryCondition right_bc,
    const T* left_tangent,
    const T* right_tangent,
    internal::StorageType<T, NumPoints*NumDims> &m
)
{
    // Initialize auxilliary storage
    internal::StorageType<T, NumPoints> a(num_points), b(num_points), c(num_points);

    for(std::size_t i = 0; i < num_points; i++)
    {
        if(i == 0)
        {
            // left boundary
            switch(left_bc)
            {
            case BoundaryCondition::Hermite:
                a[i] = 0.0;
                b[i] = 2.0;
                c[i] = 1.0;
                for(std::size_t j = 0; j < num_dims; j++)
                {
                    // store d in moments_
                    T tangent_component = 0.0;
                    if(left_tangent) tangent_component = left_tangent[j];
                    m[i*num_dims+j] = 6.0 * ((points[(i+1)*num_dims+j] - points[i*num_dims+j]) - tangent_component);
                }
                break;
            case BoundaryCondition::Periodic:
                a[i] = 1.0;
                b[i] = 4.0;
                c[i] = 1.0;
                for(std::size_t j = 0; j < num_dims; j++)
                {
                    // store d in moments_
                    m[i*num_dims+j] = 6.0 * ((points[(i+1)*num_dims+j] - points[i*num_dims+j])
                            - (points[i*num_dims+j] - points[(num_points-1)*num_dims+j]));
                }
                break;
            // TODO
            //case BoundaryCondition::NotAKnot:
            //    break;
            //
            default: // BoundaryCondition::Natural
                a[i] = 0.0;
                b[i] = 1.0;
                c[i] = 0.0;
                for(std::size_t j = 0; j < num_dims; j++)
                {
                    // store d in moments_
                    m[i*num_dims+j] = 0.0;
                }
            }
        }
        else if(i == (num_points-1))
        {
            // right boundary
            switch(right_bc)
            {
            case BoundaryCondition::Hermite:
                a[i] = 1.0;
                b[i] = 2.0;
                c[i] = 0.0;
                for(std::size_t j = 0; j < num_dims; j++)
                {
                    // store d in moments_
                    T tangent_component = 0.0;
                    if(right_tangent) tangent_component = right_tangent[j];
                    m[i*num_dims+j] = 6.0 * (tangent_component - (points[i*num_dims+j] - points[(i-1)*num_dims+j]));
                }
                break;
            case BoundaryCondition::Periodic:
                a[i] = 1.0;
                b[i] = 4.0;
                c[i] = 1.0;
                for(std::size_t j = 0; j < num_dims; j++)
                {
                    // store d in moments_
                    m[i*num_dims+j] = 6.0 * ((points[0+j] - points[(num_points-1)*num_dims+j])
                        - (points[(num_points-1)*num_dims+j] - points[(num_points-2)*num_dims+j]));
                }
                break;
            // TODO
            //case BoundaryCondition::NotAKnot:
            //    break;
            ///
            default: // BoundaryCondition::Natural
                a[i] = 0.0;
                b[i] = 1.0;
                c[i] = 0.0;
                for(std::size_t j = 0; j < num_dims; j++)
                {
                    // store d in moments_
                    m[i*num_dims+j] = 0.0;
                }
            }
        }
        else
        {
            // inner node
            a[i] = 1.0;
            b[i] = 4.0;
            c[i] = 1.0;
            for(std::size_t j = 0; j < num_dims; j++)
            {
                // store d in moments_
                m[i*num_dims+j] = 6.0 * ((points[(i+1)*num_dims+j] - points[i*num_dims+j])
                    - (points[i*num_dims+j] - points[(i-1)*num_dims+j]));
            }
        }
    }

    // Solve spline problem
    if(a[0] != 0 || c[num_points-1] != 0)
    {
        // perturbed problem
        internal::StorageType<T, NumPoints> q(num_points);
        tdma(num_points, num_dims, a, b, c, m, q.data());
    }
    else
    {
        // strictly tridiagonal problem
        tdma(num_points, num_dims, a, b, c, m);
    }
}

template<typename T, std::size_t NumPoints, std::size_t NumDims>
void Spline<T, NumPoints, NumDims>::tdma(
    const std::size_t num_points,
    const std::size_t num_dims,
    internal::StorageType<T, NumPoints> &a,
    internal::StorageType<T, NumPoints> &b,
    internal::StorageType<T, NumPoints> &c,
    internal::StorageType<T, NumPoints*NumDims> &d,
    T *q
)
{
    // Perturbed problem?
    bool is_perturbed = a[0] != 0 || c[num_points-1] != 0;
    T vn = 0.0;
    if(is_perturbed)
    {
        assert(q && "q must not be a null pointer.");

        // Initialize q and u with zero
        for(std::size_t i = 1; i <num_points; i++) q[i] = 0.0;

        // Modify problem
        vn = a[0]/b[0];
        q[0] = -b[0];
        q[num_points-1] = c[num_points-1];
        a[0] = 0;
        b[0] = 2*b[0];
        b[num_points-1] = b[num_points-1] + c[num_points-1]*vn;
        c[num_points-1] = 0;
    }

    // Forward elimination
    // i = 1 ... n:
    for(std::size_t i = 1; i < num_points; i++)
    {
        T f = a[i]/b[i-1];
        b[i] = b[i] - f*c[i-1];
        if(is_perturbed) q[i] = q[i] - f*q[i-1];
        for(std::size_t j = 0; j < num_dims; j++)
        {
            d[i*num_dims+j] = d[i*num_dims+j] - f*d[(i-1)*num_dims+j];
        }
    }

    // Backward substitution
    // i = n:
    if(is_perturbed) q[num_points-1] = q[num_points-1]/b[num_points-1];
    for(std::size_t j = 0; j < num_dims; j++)
    {
        d[(num_points-1)*num_dims+j] = d[(num_points-1)*num_dims+j]/b[num_points-1];
    }
    // i = n-1 ... 0:
    for(int i = num_points-2; i >= 0; i--)
    {
        if(is_perturbed) q[i] = (q[i] - c[i]*q[i+1])/b[i];
        for(std::size_t j = 0; j < num_dims; j++)
        {
            d[i*num_dims+j] = (d[i*num_dims+j] - c[i]*d[(i+1)*num_dims+j])/b[i];
        }
    }

    if(is_perturbed)
    {
        // Reconstruct solution
        T vq = q[0] - q[num_points-1]*vn;
        for(std::size_t j = 0; j < num_dims; j++)
        {
            T vy = d[j] - d[(num_points-1)*num_dims+j]*vn;
            T k = vy/(1 + vq);
            for(std::size_t i = 0; i < num_points; i++)
            {
                d[i*num_dims+j] = d[i*num_dims+j] - k*q[i];
            }
        }
    }
}

} // namespace: parametric_cubic_spline