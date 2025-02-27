#ifndef _MARRAY_MARRAY_VIEW_HPP_
#define _MARRAY_MARRAY_VIEW_HPP_

#include "marray_base.hpp"

namespace MArray
{

template <typename Type, int NDim>
class marray_view : public marray_base<Type, NDim, marray_view<Type, NDim>, false>
{
    protected:
        typedef marray_base<Type, NDim, marray_view, false> base;

        using base::len_;
        using base::stride_;
        using base::data_;

    public:
        using typename base::value_type;
        using typename base::pointer;
        using typename base::const_pointer;
        using typename base::reference;
        using typename base::const_reference;

        /***********************************************************************
         *
         * Constructors
         *
         **********************************************************************/

        marray_view() {}

        marray_view(const marray_view& other)
        {
            reset(other);
        }

        template <typename U, typename D, bool O,
            typename=detail::enable_if_convertible_t<
                typename marray_base<U, NDim, D, O>::cptr,pointer>>
        marray_view(const marray_base<U, NDim, D, O>& other)
        {
            reset(other);
        }

        template <typename U, typename D, bool O,
            typename=detail::enable_if_convertible_t<
                typename marray_base<U, NDim, D, O>::pointer,pointer>>
        marray_view(marray_base<U, NDim, D, O>&& other)
        {
            reset(other);
        }

        template <typename U, typename D, bool O,
            typename=detail::enable_if_convertible_t<
                typename marray_base<U, NDim, D, O>::pointer,pointer>>
        marray_view(marray_base<U, NDim, D, O>& other)
        {
            reset(other);
        }

        template <typename U, int OldNDim, int NIndexed, typename... Dims,
            typename=detail::enable_if_convertible_t<U*,pointer>>
        marray_view(const marray_slice<U, OldNDim, NIndexed, Dims...>& other)
        {
            reset(other);
        }

        marray_view(const detail::array_1d<len_type>& len, pointer ptr,
                    layout layout = DEFAULT)
        {
            reset(len, ptr, layout);
        }

        marray_view(const detail::array_1d<len_type>& len, pointer ptr,
                    const detail::array_1d<stride_type>& stride)
        {
            reset(len, ptr, stride);
        }

        /***********************************************************************
         *
         * Base operations
         *
         **********************************************************************/

        marray_view& operator=(const marray_view& other)
        {
            return base::operator=(other);
        }

        using base::operator=;
        using base::operator+=;
        using base::operator-=;
        using base::operator*=;
        using base::operator/=;
        using base::reset;
        using base::cview;
        using base::view;
        using base::cbegin;
        using base::begin;
        using base::cend;
        using base::end;
        using base::crbegin;
        using base::rbegin;
        using base::crend;
        using base::rend;
        using base::shifted;
        using base::shifted_up;
        using base::shifted_down;
        using base::permuted;
        using base::transposed;
        using base::T;
        using base::lowered;
        using base::cfront;
        using base::front;
        using base::cback;
        using base::back;
        using base::operator[];
        using base::operator();
        using base::cdata;
        using base::data;
        using base::length;
        using base::lengths;
        using base::stride;
        using base::strides;
        using base::dimension;
        using base::size;

        /***********************************************************************
         *
         * Mutating shift
         *

         **********************************************************************/

        void shift(const detail::array_1d<len_type>& n_)
        {
            MARRAY_ASSERT(n_.size() == NDim);

            std::array<len_type, NDim> n;
            n_.slurp(n);

            for (auto dim : range(NDim))
                shift(dim, n[dim]);
        }

        template <typename=void, int N=NDim, typename=detail::enable_if_t<N==1>>
        void shift(len_type n)
        {
            shift(0, n);
        }

        template <int Dim>
        void shift(len_type n)
        {
            static_assert(Dim >= 0 && Dim < NDim, "Dim out of range");
            shift(Dim, n);
        }

        void shift(int dim, len_type n)
        {
            MARRAY_ASSERT(dim >= 0 && dim < NDim);
            data_ += n*stride_[dim];
        }

        template <typename=void, int N=NDim, typename=detail::enable_if_t<N==1>>
        void shift_down()
        {
            shift_down(0);
        }

        template <int Dim>
        void shift_down()
        {
            shift_down(Dim);
        }

        void shift_down(int dim)
        {
            shift(dim, len_[dim]);
        }

        template <int Dim>
        void shift_up()
        {
            shift_up(Dim);
        }

        template <typename=void, int N=NDim, typename=detail::enable_if_t<N==1>>
        void shift_up()
        {
            shift_up(0);
        }

        void shift_up(int dim)
        {
            shift(dim, -len_[dim]);
        }

        /***********************************************************************
         *
         * Mutating permute
         *
         **********************************************************************/

        void permute(const detail::array_1d<int>& perm_)
        {
            MARRAY_ASSERT(perm_.size() == NDim);

            std::array<len_type, NDim> len = len_;
            std::array<stride_type, NDim> stride = stride_;
            std::array<int, NDim> perm;
            perm_.slurp(perm);

            for (auto i : range(NDim))
            {
                MARRAY_ASSERT(perm[i] < NDim);
                for (auto j : range(i+1,NDim))
                    MARRAY_ASSERT(perm[i] != perm[j]);

                len_[i] = len[perm[i]];
                stride_[i] = stride[perm[i]];
            }
        }

        template <int N=NDim, typename=detail::enable_if_t<N==2>>
        void transpose()
        {
            permute({1, 0});
        }

        /***********************************************************************
         *
         * Mutating reversal
         *
         **********************************************************************/

        void reverse()
        {
            for (auto i : range(NDim)) reverse(i);
        }

        template <int Dim>
        void reverse()
        {
            reverse(Dim);
        }

        void reverse(int dim)
        {
            MARRAY_ASSERT(dim >= 0 && dim < NDim);
            data_ += (len_[dim]-1)*stride_[dim];
            stride_[dim] = -stride_[dim];
        }

        /***********************************************************************
         *
         * Basic setters
         *
         **********************************************************************/

        pointer data(pointer ptr)
        {
            std::swap(ptr, data_);
            return ptr;
        }

        template <int Dim>
        len_type length(len_type len)
        {
            static_assert(Dim >= 0 && Dim < NDim, "Dim out of range");
            return length(Dim, len);
        }

        len_type length(int dim, len_type len)
        {
            MARRAY_ASSERT(dim >= 0 && dim < NDim);
            std::swap(len, len_[dim]);
            return len;
        }

        template <int Dim>
        stride_type stride(stride_type s)
        {
            static_assert(Dim >= 0 && Dim < NDim, "Dim out of range");
            return stride(Dim, s);
        }

        stride_type stride(int dim, stride_type s)
        {
            MARRAY_ASSERT(dim >= 0 && dim < NDim);
            std::swap(s, stride_[dim]);
            return s;
        }

        /***********************************************************************
         *
         * Swap
         *
         **********************************************************************/

        void swap(marray_view& other)
        {
            base::swap(other);
        }

        friend void swap(marray_view& a, marray_view& b)
        {
            a.swap(b);
        }
};

}

#endif
