#ifndef _MARRAY_MARRAY_BASE_HPP_
#define _MARRAY_MARRAY_BASE_HPP_

#include "utility.hpp"
#include "range.hpp"

namespace MArray
{

template <typename Type, int NDim, int NIndexed, typename... Dims>
class marray_slice;

template <typename Type, int NDim, typename Derived, bool Owner>
class marray_base;

template <typename Type, int NDim>
class marray_view;

template <typename Type, int NDim, typename Allocator=std::allocator<Type>>
class marray;

template <typename Type, typename Derived, bool Owner>
class varray_base;

template <typename Type>
class varray_view;

template <typename Type, typename Allocator=std::allocator<Type>>
class varray;

template <typename Expr>
struct is_expression_arg_or_scalar;

/*
 * Convenient names for 1- and 2-dimensional array types.
 */
template <typename Type> using row_view = marray_view<Type, 1>;
template <typename Type, typename Allocator=std::allocator<Type>> using row = marray<Type, 1, Allocator>;

template <typename Type> using matrix_view = marray_view<Type, 2>;
template <typename Type, typename Allocator=std::allocator<Type>> using matrix = marray<Type, 2, Allocator>;

namespace detail
{

struct len_type_wrapper
{
    len_type value;

    len_type_wrapper(  signed     short value) : value(value) {}
    len_type_wrapper(unsigned     short value) : value(value) {}
    len_type_wrapper(  signed       int value) : value(value) {}
    len_type_wrapper(unsigned       int value) : value(value) {}
    len_type_wrapper(  signed      long value) : value(value) {}
    len_type_wrapper(unsigned      long value) : value(value) {}
    len_type_wrapper(  signed long long value) : value(value) {}
    len_type_wrapper(unsigned long long value) : value(value) {}

    operator len_type() const { return value; }
};

typedef std::initializer_list<len_type_wrapper> len_type_init;

template <typename Type, int NDim>
struct initializer_type;

template <typename Type>
struct initializer_type<Type, 0u>
{
    typedef Type type;
};

template <typename Type, int NDim>
struct initializer_type
{
    typedef std::initializer_list<
        typename initializer_type<Type, NDim-1>::type> type;
};

template <typename Func, typename Arg, typename... Args>
detail::enable_if_t<sizeof...(Args) &&
                    std::is_same<decltype(std::declval<Func&&>()(
                        std::declval<Arg&&>(),
                        std::declval<Args&&>()...)),void>::value>
call(Func&& f, Arg&& arg, Args&&... args)
{
    f(std::forward<Arg>(arg), std::forward<Args>(args)...);
}

template <typename Func, typename Arg, typename... Args>
detail::enable_if_t<std::is_same<decltype(std::declval<Func&&>()(
                        std::declval<Arg&&>())),void>::value>
call(Func&& f, Arg&& arg, Args&&...)
{
    f(std::forward<Arg>(arg));
}

template <typename Array>
class marray_iterator
{
    protected:
        Array* array_ = nullptr;
        len_type i_ = 0;

    public:
        typedef std::random_access_iterator_tag iterator_category;
        typedef decltype((*array_)[0]) value_type;
        typedef len_type difference_type;
        typedef typename std::remove_reference<value_type>::type* pointer;
        typedef value_type& reference;

        marray_iterator(Array& array, len_type i)
        : array_(&array), i_(i) {}

        bool operator==(const marray_iterator& other) const
        {
            return i_ == other.i_;
        }

        bool operator!=(const marray_iterator& other) const
        {
            return !(*this == other);
        }

        value_type operator*() const
        {
            return (*array_)[i_];
        }

        marray_iterator& operator++()
        {
            i_++;
            return *this;
        }

        marray_iterator operator++(int)
        {
            marray_iterator old(*this);
            i_++;
            return old;
        }

        marray_iterator& operator--()
        {
            i_--;
            return *this;
        }

        marray_iterator operator--(int)
        {
            marray_iterator old(*this);
            i_--;
            return old;
        }

        marray_iterator& operator+=(difference_type n)
        {
            i_ += n;
            return *this;
        }

        marray_iterator operator+(difference_type n) const
        {
            return marray_iterator(*array_, i_+n);
        }

        friend marray_iterator operator+(difference_type n, const marray_iterator& i)
        {
            return marray_iterator(*i.array_, i.i_+n);
        }

        marray_iterator& operator-=(difference_type n)
        {
            i_ -= n;
            return *this;
        }

        marray_iterator operator-(difference_type n) const
        {
            return marray_iterator(*array_, i_-n);
        }

        difference_type operator-(const marray_iterator& other) const
        {
            return i_ - other.i_;
        }

        bool operator<(const marray_iterator& other) const
        {
            return i_ < other.i_;
        }

        bool operator<=(const marray_iterator& other) const
        {
            return !(other < *this);
        }

        bool operator>(const marray_iterator& other) const
        {
            return other < *this;
        }

        bool operator>=(const marray_iterator& other) const
        {
            return !(*this < other);
        }

        value_type operator[](difference_type n) const
        {
            return (*array_)[i_+n];
        }

        friend void swap(marray_iterator& a, marray_iterator& b)
        {
            using std::swap;
            swap(a.array_, b.array_);
            swap(a.i_, b.i_);
        }
};

template <typename T>
struct is_row : std::false_type {};

template <typename T, typename D, bool O>
struct is_row<marray_base<T, 1, D, O>> : std::true_type {};

template <typename T>
struct is_row<marray_view<T, 1>> : std::true_type {};

template <typename T, typename A>
struct is_row<marray<T, 1, A>> : std::true_type {};

template <typename T, typename U, typename=void>
struct is_row_of : std::false_type {};

template <typename T, typename U>
struct is_row_of<T, U, enable_if_t<is_row<T>::value &&
    std::is_assignable<U&,typename T::value_type>::value>> : std::true_type {};

template <typename T>
struct is_matrix : std::false_type {};

template <typename T, typename D, bool O>
struct is_matrix<marray_base<T, 2, D, O>> : std::true_type {};

template <typename T>
struct is_matrix<marray_view<T, 2>> : std::true_type {};

template <typename T, typename A>
struct is_matrix<marray<T, 2, A>> : std::true_type {};

template <typename T, typename U, typename=void>
struct is_matrix_of : std::false_type {};

template <typename T, typename U>
struct is_matrix_of<T, U, enable_if_t<is_matrix<T>::value &&
    std::is_assignable<U&,typename T::value_type>::value>> : std::true_type {};

template <typename T, typename U, typename V=void>
using enable_if_matrix_of_t = enable_if_t<is_matrix_of<T,U>::value, V>;

template <typename T, typename U>
struct is_1d_container_of :
    std::integral_constant<bool, is_row_of<T,U>::value ||
                                 is_container_of<T,U>::value> {};

template <typename T, typename U, typename V=void>
using enable_if_1d_container_of_t = enable_if_t<is_1d_container_of<T,U>::value, V>;

template <typename T, typename U>
struct is_2d_container_of :
    std::integral_constant<bool, is_matrix_of<T,U>::value ||
                                 is_container_of_containers_of<T,U>::value> {};

template <typename T, typename U, typename V=void>
using enable_if_2d_container_of_t = enable_if_t<is_2d_container_of<T,U>::value, V>;

template <typename T>
enable_if_t<is_container<T>::value,len_type>
length(const T& len)
{
    return len.size();
}

template <typename T>
enable_if_t<is_row<T>::value,len_type>
length(const T& len)
{
    return len.length(0);
}

template <typename T>
enable_if_t<is_container_of_containers<T>::value,len_type>
length(const T& len, int dim)
{
    if (dim == 0) return len.size();
    else
    {
        MARRAY_ASSERT(dim == 1);
        auto it = len.begin();
        if (it == len.end()) return 0;
        len_type l = it->size();
        while (++it != len.end()) MARRAY_ASSERT((len_type)it->size() == l);
        return l;
    }
}

template <typename T>
enable_if_t<is_matrix<T>::value,len_type>
length(const T& len, int dim)
{
    return len.length(dim);
}

template <typename T>
class array_1d
{
    protected:
        struct adaptor_base
        {
            len_type len;

            adaptor_base(len_type len) : len(len) {}

            virtual ~adaptor_base() {}

            virtual void slurp(T*) const = 0;

            virtual adaptor_base& move(adaptor_base& other) = 0;
        };

        template <typename U>
        struct adaptor : adaptor_base
        {
            U data;

            adaptor(U data)
            : adaptor_base(detail::length(data)), data(std::move(data)) {}

            virtual void slurp(T* x) const override
            {
                std::copy_n(data.begin(), this->len, x);
            }

            virtual adaptor_base& move(adaptor_base& other) override
            {
            	return *(new (static_cast<adaptor*>(&other)) adaptor(std::move(*this)));
            }
        };

        constexpr static size_t _s1 = sizeof(adaptor<T*&>);
        constexpr static size_t _s2 = sizeof(adaptor<short_vector<T,MARRAY_OPT_NDIM>>);
        constexpr static size_t max_adaptor_size = _s1 > _s2 ? _s1 : _s2;

        std::aligned_storage_t<max_adaptor_size> raw_adaptor_;
        adaptor_base& adaptor_;

        template <typename U>
        adaptor_base& adapt(U&& data)
        {
            return *(new (&raw_adaptor_) adaptor<U>(data));
        }

    public:
        array_1d()
        : adaptor_(adapt(std::array<T,0>{})) {}

        array_1d(const array_1d& other)
        : adaptor_(other.adaptor_.move(reinterpret_cast<adaptor_base&>(raw_adaptor_))) {}

        template <typename... Args, typename =
            detail::enable_if_t<detail::are_convertible<T,Args...>::value>>
        array_1d(Args&&... args)
        : adaptor_(adapt(short_vector<T,MARRAY_OPT_NDIM>{(T)args...})) {}

        template <typename U, typename=detail::enable_if_1d_container_of_t<U,T>>
        array_1d(const U& data)
        : adaptor_(adapt(data)) {}

        ~array_1d() { adaptor_.~adaptor_base(); }

        template <size_t N>
        void slurp(std::array<T, N>& x) const
        {
            MARRAY_ASSERT((len_type)N >= size());
            adaptor_.slurp(x.data());
        }

        void slurp(std::vector<T>& x) const
        {
            x.resize(size());
            adaptor_.slurp(x.data());
        }

        template <size_t N>
        void slurp(short_vector<T, N>& x) const
        {
            x.resize(size());
            adaptor_.slurp(x.data());
        }

        void slurp(row<T>& x) const
        {
            x.reset({size()});
            adaptor_.slurp(x.data());
        }

        len_type size() const { return adaptor_.len; }
};

template <typename T>
class array_2d
{
    protected:
        struct adaptor_base
        {
            std::array<len_type,2> len;

            adaptor_base(len_type len0, len_type len1) : len{len0, len1} {}

            virtual ~adaptor_base() {}

            virtual void slurp(T* x, len_type rs, len_type cs) const = 0;

            virtual void slurp(std::vector<std::vector<T>>&) const = 0;
        };

        template <typename U>
        struct adaptor : adaptor_base
        {
            static constexpr bool IsMatrix = is_matrix<typename std::decay<U>::type>::value;

            U data;

            adaptor(U data)
            : adaptor_base(detail::length(data, 0), detail::length(data, 1)),
              data(data) {}

            template <bool IsMatrix_ = IsMatrix>
            typename std::enable_if<!IsMatrix_>::type
            do_slurp(T* x, len_type rs, len_type cs) const
            {
                int i = 0;
                for (auto it = this->data.begin(), end = this->data.end();it != end;++it)
                {
                    int j = 0;
                    for (auto it2 = it->begin(), end2 = it->end();it2 != end2;++it2)
                    {
                        x[i*rs + j*cs] = *it2;
                        j++;
                    }
                    i++;
                }
            }

            template <bool IsMatrix_ = IsMatrix>
            typename std::enable_if<!IsMatrix_>::type
            do_slurp(std::vector<std::vector<T>>& x) const
            {
                x.clear();
                for (auto it = this->data.begin(), end = this->data.end();it != end;++it)
                {
                    x.emplace_back(it->begin(), it->end());
                }
            }

            template <bool IsMatrix_ = IsMatrix>
            typename std::enable_if<IsMatrix_>::type
            do_slurp(T* x, len_type rs, len_type cs) const
            {
                for (len_type i = 0;i < this->len[0];i++)
                {
                    for (len_type j = 0;j < this->len[1];j++)
                    {
                        x[i*rs + j*cs] = this->data[i][j];
                    }
                }
            }

            template <bool IsMatrix_ = IsMatrix>
            typename std::enable_if<IsMatrix_>::type
            do_slurp(std::vector<std::vector<T>>& x) const
            {
                x.resize(this->len[0]);
                for (len_type i = 0;i < this->len[0];i++)
                {
                    x[i].resize(this->len[1]);
                    for (len_type j = 0;j < this->len[1];j++)
                    {
                        x[i][j] = this->data[i][j];
                    }
                }
            }

            virtual void slurp(T* x, len_type rs, len_type cs) const override
            {
                do_slurp(x, rs, cs);
            }

            virtual void slurp(std::vector<std::vector<T>>& x) const override
            {
                do_slurp(x);
            }
        };

        constexpr static size_t _s1 = sizeof(adaptor<std::initializer_list<std::initializer_list<int>>>);
        constexpr static size_t _s2 = sizeof(adaptor<const matrix<int>&>);
        constexpr static size_t max_adaptor_size = _s1 > _s2 ? _s1 : _s2;

        std::aligned_storage_t<max_adaptor_size> raw_adaptor_;
        adaptor_base& adaptor_;

        template <typename U>
        adaptor_base& adapt(U data)
        {
            return *(new (&raw_adaptor_) adaptor<U>(data));
        }

        adaptor_base& adapt(const array_2d& other)
        {
            memcpy(&raw_adaptor_, &other.raw_adaptor_, sizeof(raw_adaptor_));

            return *reinterpret_cast<adaptor_base*>(
                reinterpret_cast<char*>(this) +
                (reinterpret_cast<const char*>(&other.adaptor_) -
                 reinterpret_cast<const char*>(&other)));
        }

    public:
        array_2d(const array_2d& other)
        : adaptor_(adapt(other)) {}

        array_2d(std::initializer_list<std::initializer_list<T>> data)
        : adaptor_(adapt<std::initializer_list<std::initializer_list<T>>>(data)) {}

        template <typename U, typename=detail::enable_if_assignable_t<T&,U>>
        array_2d(std::initializer_list<std::initializer_list<U>> data)
        : adaptor_(adapt<std::initializer_list<std::initializer_list<U>>>(data)) {}

        template <typename U, typename=detail::enable_if_1d_container_of_t<U,T>>
        array_2d(std::initializer_list<U> data)
        : adaptor_(adapt<std::initializer_list<U>>(data)) {}

        template <typename U, typename=detail::enable_if_2d_container_of_t<U,T>>
        array_2d(const U& data)
        : adaptor_(adapt<const U&>(data)) {}

        ~array_2d() { adaptor_.~adaptor_base(); }

        void slurp(std::vector<std::vector<T>>& x) const { adaptor_.slurp(x); }

        template <size_t M, size_t N>
        void slurp(std::array<std::array<T, N>, M>& x) const
        {
            MARRAY_ASSERT((len_type)M >= length(0));
            MARRAY_ASSERT((len_type)N >= length(1));
            adaptor_.slurp(&x[0][0], N, 1);
        }

        template <size_t M, size_t N>
        void slurp(short_vector<std::array<T, N>, M>& x) const
        {
            x.resize(length(0));
            MARRAY_ASSERT((len_type)N >= length(1));
            adaptor_.slurp(&x[0][0], N, 1);
        }

        void slurp(matrix<T>& x, layout layout = DEFAULT) const
        {
            x.reset({length(0), length(1)}, layout);
            adaptor_.slurp(x.data(), x.stride(0), x.stride(1));
        }

        len_type length(int dim) const
        {
            MARRAY_ASSERT(dim >= 0 && dim < 2);
            return adaptor_.len[dim];
        }
};

}
}

#include "expression.hpp"
#include "marray_slice.hpp"
#include "miterator.hpp"

namespace MArray
{

template <typename Type, int NDim, typename Derived, bool Owner>
class marray_base
{
    static_assert(NDim > 0, "NDim must be positive");

    template <typename, int, typename, bool> friend class marray_base;
    template <typename, int> friend class marray_view;
    template <typename, int, typename> friend class marray;
    template <typename, typename, bool> friend class varray_base;
    template <typename> friend class varray_view;
    template <typename, typename> friend class varray;

    public:
        typedef Type value_type;
        typedef Type* pointer;
        typedef const Type* const_pointer;
        typedef Type& reference;
        typedef const Type& const_reference;
        typedef typename detail::initializer_type<Type, NDim>::type
            initializer_type;
        typedef detail::marray_iterator<marray_base> iterator;
        typedef detail::marray_iterator<const marray_base> const_iterator;
        typedef std::reverse_iterator<iterator> reverse_iterator;
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

        typedef typename std::conditional<Owner,const Type,Type>::type ctype;
        typedef ctype& cref;
        typedef ctype* cptr;

    protected:
        std::array<len_type, NDim> len_ = {};
        std::array<stride_type, NDim> stride_ = {};
        pointer data_ = nullptr;

        /***********************************************************************
         *
         * Reset
         *
         **********************************************************************/

        void reset()
        {
            data_ = nullptr;
            len_ = {};
            stride_ = {};
        }

        template <typename U, bool O, typename D,
            typename=detail::enable_if_convertible_t<
                typename marray_base<U, NDim, D, O>::cptr,pointer>>
        void reset(const marray_base<U, NDim, D, O>& other)
        {
            reset(const_cast<marray_base<U, NDim, D, O>&>(other));
        }

        template <typename U, bool O, typename D,
            typename=detail::enable_if_convertible_t<
                typename marray_base<U, NDim, D, O>::pointer,pointer>>
        void reset(marray_base<U, NDim, D, O>& other)
        {
            data_ = other.data_;
            len_ = other.len_;
            stride_ = other.stride_;
        }

        template <typename U, int OldNDim, int NIndexed,
            typename... Dims, typename=detail::enable_if_convertible_t<U*,pointer>>
        void reset(const marray_slice<U, OldNDim, NIndexed, Dims...>& other)
        {
            reset(other.view());
        }

        void reset(const detail::array_1d<len_type>& len, pointer ptr,
                   layout layout = DEFAULT)
        {
            reset(len, ptr, strides(len, layout));
        }

        void reset(const detail::array_1d<len_type>& len, pointer ptr,
                   const detail::array_1d<stride_type>& stride)
        {
            data_ = ptr;
            len.slurp(len_);
            stride.slurp(stride_);
        }

        void reset(initializer_type data, layout layout = DEFAULT)
        {
            set_lengths(0, len_, data);
            stride_ = strides(len_, layout);
            set_data(0, data_, data);
        }

        /***********************************************************************
         *
         * Private helper functions
         *
         **********************************************************************/

        template <typename Ptr, typename Func, int... I>
        void for_each_element(Func&& f,
                              detail::integer_sequence<int, I...>) const
        {
            miterator<NDim, 1> it(len_, stride_);
            Ptr ptr = const_cast<Ptr>(data_);
            while (it.next(ptr))
                detail::call(std::forward<Func>(f), *ptr, it.position()[I]...);
        }

        void set_lengths(int i, std::array<len_type, NDim>& len,
                         std::initializer_list<Type> data)
        {
            len[i] = data.size();
        }

        template <typename U>
        void set_lengths(int i, std::array<len_type, NDim>& len,
                         std::initializer_list<U> data)
        {
            len[i] = data.size();
            set_lengths(i+1, len, *data.begin());
        }

        void set_data(int i, pointer ptr, std::initializer_list<Type> data)
        {
            auto it = data.begin();
            for (len_type j = 0;j < len_[i];j++)
            {
                ptr[j*stride_[i]] = *it;
                ++it;
            }
        }

        template <typename U>
        void set_data(int i, pointer ptr, std::initializer_list<U> data)
        {
            auto it = data.begin();
            for (len_type j = 0;j < len_[i];j++)
            {
                set_data(i+1, ptr + j*stride_[i], *it);
                ++it;
            }
        }

        void swap(marray_base& other)
        {
            using std::swap;
            swap(data_, other.data_);
            swap(len_, other.len_);
            swap(stride_, other.stride_);
        }

    public:

        /***********************************************************************
         *
         * Static helper functions
         *
         **********************************************************************/


        static std::array<stride_type, NDim>
        strides(const detail::array_1d<len_type>& len_, layout layout = DEFAULT)
        {
            //TODO: add alignment option

            MARRAY_ASSERT(len_.size() == NDim);

            std::array<len_type, NDim> len;
            len_.slurp(len);
            std::array<stride_type, NDim> stride;

            if (layout == ROW_MAJOR)
            {
                stride[NDim-1] = 1;
                for (auto i : reversed_range(NDim-1))
                    stride[i] = stride[i+1]*len[i+1];
            }
            else
            {
                stride[0] = 1;
                for (auto i : range(1,NDim))
                    stride[i] = stride[i-1]*len[i-1];
            }

            return stride;
        }

        static stride_type size(const detail::array_1d<len_type>& len_)
        {
            //TODO: add alignment option

            len_vector len;
            len_.slurp(len);

            stride_type s = 1;
            for (auto& l : len) s *= l;
            return s;
        }

        /***********************************************************************
         *
         * Operators
         *
         **********************************************************************/

        Derived& operator=(const marray_base& other)
        {
            return operator=<>(other);
        }

        Derived& operator=(initializer_type data)
        {
            std::array<len_type, NDim> len;
            set_lengths(0, len, data);
            MARRAY_ASSERT(len == len_);
            set_data(0, data_, data);
            return static_cast<Derived&>(*this);
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        Derived& operator=(const Expression& other)
        {
            assign_expr(*this, other);
            return static_cast<Derived&>(*this);
        }

        template <typename Expression, bool O=Owner,
            typename=detail::enable_if_t<!O && is_expression_arg_or_scalar<Expression>::value>>
        const Derived& operator=(const Expression& other) const
        {
            assign_expr(*this, other);
            return static_cast<const Derived&>(*this);
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        Derived& operator+=(const Expression& other)
        {
            *this = *this + other;
            return static_cast<Derived&>(*this);
        }

        template <typename Expression, bool O=Owner,
            typename=detail::enable_if_t<!O && is_expression_arg_or_scalar<Expression>::value>>
        const Derived& operator+=(const Expression& other) const
        {
            *this = *this + other;
            return static_cast<const Derived&>(*this);
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        Derived& operator-=(const Expression& other)
        {
            *this = *this - other;
            return static_cast<Derived&>(*this);
        }

        template <typename Expression, bool O=Owner,
            typename=detail::enable_if_t<!O && is_expression_arg_or_scalar<Expression>::value>>
        const Derived& operator-=(const Expression& other) const
        {
            *this = *this - other;
            return static_cast<const Derived&>(*this);
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        Derived& operator*=(const Expression& other)
        {
            *this = *this * other;
            return static_cast<Derived&>(*this);
        }

        template <typename Expression, bool O=Owner,
            typename=detail::enable_if_t<!O && is_expression_arg_or_scalar<Expression>::value>>
        const Derived& operator*=(const Expression& other) const
        {
            *this = *this * other;
            return static_cast<const Derived&>(*this);
        }

        template <typename Expression,
            typename=detail::enable_if_t<is_expression_arg_or_scalar<Expression>::value>>
        Derived& operator/=(const Expression& other)
        {
            *this = *this / other;
            return static_cast<Derived&>(*this);
        }

        template <typename Expression, bool O=Owner,
            typename=detail::enable_if_t<!O && is_expression_arg_or_scalar<Expression>::value>>
        const Derived& operator/=(const Expression& other) const
        {
            *this = *this / other;
            return static_cast<const Derived&>(*this);
        }

        template <typename U, int N, typename D, bool O>
        detail::enable_if_t<N==NDim, bool>
        operator==(const marray_base<U, N, D, O>& other) const
        {
            if (len_ != other.len_) return false;

            miterator<NDim, 2> it(len_, stride_, other.stride_);

            auto a = data_;
            auto b = other.data_;
            while (it.next(a, b))
            {
                if (*a != *b) return false;
            }

            return true;
        }

        template <typename U, int N, typename D, bool O>
        detail::enable_if_t<N!=NDim, bool>
        operator==(const marray_base<U, N, D, O>&) const
        {
            return false;
        }

        template <typename U, int N, typename D, bool O>
        bool operator!=(const marray_base<U, N, D, O>& other) const
        {
            return !(*this == other);
        }

        /***********************************************************************
         *
         * Views
         *
         **********************************************************************/

        marray_view<const Type, NDim> cview() const
        {
            return const_cast<marray_base&>(*this).view();
        }

        marray_view<ctype, NDim> view() const
        {
            return const_cast<marray_base&>(*this).view();
        }

        marray_view<Type, NDim> view()
        {
            return *this;
        }

        friend marray_view<const Type, NDim> cview(const marray_base& x)
        {
            return x.view();
        }

        friend marray_view<ctype, NDim> view(const marray_base& x)
        {
            return x.view();
        }

        friend marray_view<Type, NDim> view(marray_base& x)
        {
            return x.view();
        }

        /***********************************************************************
         *
         * Iterators
         *
         **********************************************************************/

        const_iterator cbegin() const
        {
            return const_iterator{*this, 0};
        }

        const_iterator begin() const
        {
            return const_iterator{*this, 0};
        }

        iterator begin()
        {
            return iterator{*this, 0};
        }

        const_iterator cend() const
        {
            return const_iterator{*this, len_[0]};
        }

        const_iterator end() const
        {
            return const_iterator{*this, len_[0]};
        }

        iterator end()
        {
            return iterator{*this, len_[0]};
        }

        const_reverse_iterator crbegin() const
        {
            return const_reverse_iterator{end()};
        }

        const_reverse_iterator rbegin() const
        {
            return const_reverse_iterator{end()};
        }

        reverse_iterator rbegin()
        {
            return reverse_iterator{end()};
        }

        const_reverse_iterator crend() const
        {
            return const_reverse_iterator{begin()};
        }

        const_reverse_iterator rend() const
        {
            return const_reverse_iterator{begin()};
        }

        reverse_iterator rend()
        {
            return reverse_iterator{begin()};
        }

        /***********************************************************************
         *
         * Shifting
         *
         **********************************************************************/

        marray_view<ctype, NDim> shifted(const detail::array_1d<len_type>& n) const
        {
            return const_cast<marray_base&>(*this).shifted(n);
        }

        marray_view<Type, NDim> shifted(const detail::array_1d<len_type>& n)
        {
            marray_view<Type,NDim> r(*this);
            r.shift(n);
            return r;
        }

        template <typename=void, int N=NDim, typename=detail::enable_if_t<N==1>>
        marray_view<ctype,1> shifted(len_type n) const
        {
            return const_cast<marray_base&>(*this).shifted(n);
        }

        template <typename=void, int N=NDim, typename=detail::enable_if_t<N==1>>
        marray_view<Type,1> shifted(len_type n)
        {
            return shifted(0, n);
        }

        template <int Dim>
        marray_view<ctype, NDim> shifted(len_type n) const
        {
            return const_cast<marray_base&>(*this).shifted<Dim>(n);
        }

        template <int Dim>
        marray_view<Type, NDim> shifted(len_type n)
        {
            return shifted(Dim, n);
        }

        marray_view<ctype, NDim> shifted(int dim, len_type n) const
        {
            return const_cast<marray_base&>(*this).shifted(dim, n);
        }

        marray_view<Type, NDim> shifted(int dim, len_type n)
        {
            marray_view<Type,NDim> r(*this);
            r.shift(dim, n);
            return r;
        }

        template <typename=void, int N=NDim, typename=detail::enable_if_t<N==1>>
        marray_view<ctype,1> shifted_down() const
        {
            return const_cast<marray_base&>(*this).shifted_down();
        }

        template <typename=void, int N=NDim, typename=detail::enable_if_t<N==1>>
        marray_view<Type,1> shifted_down()
        {
            return shifted_down(0);
        }

        template <int Dim>
        marray_view<ctype,NDim> shifted_down() const
        {
            return const_cast<marray_base&>(*this).shifted_down<Dim>();
        }

        template <int Dim>
        marray_view<Type,NDim> shifted_down()
        {
            return shifted_down(Dim);
        }

        marray_view<ctype,NDim> shifted_down(int dim) const
        {
            return const_cast<marray_base&>(*this).shifted_down(dim);
        }

        marray_view<Type,NDim> shifted_down(int dim)
        {
            return shifted(dim, len_[dim]);
        }

        template <typename=void, int N=NDim, typename=detail::enable_if_t<N==1>>
        marray_view<ctype,1> shifted_up() const
        {
            return const_cast<marray_base&>(*this).shifted_up();
        }

        template <typename=void, int N=NDim, typename=detail::enable_if_t<N==1>>
        marray_view<Type,1> shifted_up()
        {
            return shifted_up(0);
        }

        template <int Dim>
        marray_view<ctype,NDim> shifted_up() const
        {
            return const_cast<marray_base&>(*this).shifted_up<Dim>();
        }

        template <int Dim>
        marray_view<Type,NDim> shifted_up()
        {
            return shifted_up(Dim);
        }

        marray_view<ctype,NDim> shifted_up(int dim) const
        {
            return const_cast<marray_base&>(*this).shifted_up(dim);
        }

        marray_view<Type,NDim> shifted_up(int dim)
        {
            return shifted(dim, -len_[dim]);
        }

        /***********************************************************************
         *
         * Permutation
         *
         **********************************************************************/

        marray_view<ctype,NDim> permuted(const detail::array_1d<int>& perm) const
        {
            return const_cast<marray_base&>(*this).permuted(perm);
        }

        marray_view<Type,NDim> permuted(const detail::array_1d<int>& perm)
        {
            marray_view<Type,NDim> r(*this);
            r.permute(perm);
            return r;
        }

        template <int N=NDim, typename=detail::enable_if_t<N==2>>
        marray_view<ctype, NDim> transposed() const
        {
            return const_cast<marray_base&>(*this).transposed();
        }

        template <int N=NDim, typename=detail::enable_if_t<N==2>>
        marray_view<Type, NDim> transposed()
        {
            return permuted({1, 0});
        }

        template <int N=NDim, typename=detail::enable_if_t<N==2>>
        marray_view<ctype, NDim> T() const
        {
            return const_cast<marray_base&>(*this).T();
        }

        template <int N=NDim, typename=detail::enable_if_t<N==2>>
        marray_view<Type, NDim> T()
        {
            return transposed();
        }

        /***********************************************************************
         *
         * Dimension change
         *
         **********************************************************************/

        template <int NewNDim>
        marray_view<ctype, NewNDim> lowered(const detail::array_1d<int>& split) const
        {
            return const_cast<marray_base&>(*this).lowered<NewNDim>(split);
        }

        template <int NewNDim>
        marray_view<Type, NewNDim> lowered(const detail::array_1d<int>& split_)
        {
            static_assert(NewNDim > 0 && NewNDim <= NDim,
                          "Cannot split into this number of dimensions");

            constexpr auto NSplit = NewNDim-1;
            MARRAY_ASSERT(split_.size() == NSplit);

            std::array<int, NSplit> split;
            split_.slurp(split);

            for (auto i : range(NSplit))
            {
                MARRAY_ASSERT(split[i] <= NDim);
                if (i != 0) MARRAY_ASSERT(split[i-1] <= split[i]);
            }

            std::array<len_type, NSplit+1> newlen;
            std::array<stride_type, NSplit+1> newstride;

            for (auto i : range(NSplit+1))
            {
                auto begin = (i == 0 ? 0 : split[i-1]);
                auto end = (i == NSplit ? NDim-1 : split[i]-1);
                if (begin > end) continue;

                if (stride_[begin] < stride_[end] ||
                    (stride_[begin] == stride_[end] && len_[begin] == 1))
                {
                    newlen[i] = len_[end];
                    newstride[i] = stride_[begin];
                    for (auto j : range(begin,end))
                    {
                        MARRAY_ASSERT(stride_[j+1] == stride_[j]*len_[j]);
                        newlen[i] *= len_[j];
                    }
                }
                else
                {
                    newlen[i] = len_[end];
                    newstride[i] = stride_[end];
                    for (auto j : range(begin,end))
                    {
                        MARRAY_ASSERT(stride_[j] == stride_[j+1]*len_[j+1]);
                        newlen[i] *= len_[j];
                    }
                }
            }

            return {newlen, data_, newstride};
        }

        /***********************************************************************
         *
         * Reversal
         *
         **********************************************************************/

        marray_view<ctype, NDim> reversed() const
        {
            return const_cast<marray_base&>(*this).reversed();
        }

        marray_view<Type, NDim> reversed()
        {
            marray_view<Type,NDim> r(*this);
            r.reverse();
            return r;
        }

        template <int Dim>
        marray_view<ctype, NDim> reversed() const
        {
            return const_cast<marray_base&>(*this).reversed<Dim>();
        }

        template <int Dim>
        marray_view<Type, NDim> reversed()
        {
            return reversed(Dim);
        }

        marray_view<ctype, NDim> reversed(int dim) const
        {
            return const_cast<marray_base&>(*this).reversed(dim);
        }

        marray_view<Type, NDim> reversed(int dim)
        {
            marray_view<Type,NDim> r(*this);
            r.reverse(dim);
            return r;
        }

        /***********************************************************************
         *
         * Slices
         *
         **********************************************************************/

        template <int Dim=0, int N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cfront() const
        {
            return const_cast<marray_base&>(*this).front<Dim>();
        }

        template <int Dim=0, int N=NDim>
        detail::enable_if_t<N==1, cref>
        front() const
        {
            return const_cast<marray_base&>(*this).front<Dim>();
        }

        template <int Dim=0, int N=NDim>
        detail::enable_if_t<N==1, reference>
        front()
        {
            static_assert(Dim == 0, "Dim out of range");
            MARRAY_ASSERT(len_[0] > 0);
            return data_[0];
        }

        template <int N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cfront(int dim) const
        {
            return const_cast<marray_base&>(*this).front(dim);
        }

        template <int N=NDim>
        detail::enable_if_t<N==1, cref>
        front(int dim) const
        {
            return const_cast<marray_base&>(*this).front(dim);
        }

        template <int N=NDim>
        detail::enable_if_t<N==1, reference>
        front(int dim)
        {
            (void)dim;
            MARRAY_ASSERT(dim == 0);
            return front();
        }

        template <int Dim, int N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cfront() const
        {
            return const_cast<marray_base&>(*this).front<Dim>();
        }

        template <int Dim, int N=NDim>
        detail::enable_if_t<N!=1, marray_view<ctype, NDim-1>>
        front() const
        {
            return const_cast<marray_base&>(*this).front<Dim>();
        }

        template <int Dim, int N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        front()
        {
            return front<N>(Dim);
        }

        template <int N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cfront(int dim) const
        {
            return const_cast<marray_base&>(*this).front(dim);
        }

        template <int N=NDim>
        detail::enable_if_t<N!=1, marray_view<ctype, NDim-1>>
        front(int dim) const
        {
            return const_cast<marray_base&>(*this).front(dim);
        }

        template <int N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        front(int dim)
        {
            MARRAY_ASSERT(dim >= 0 && dim < NDim);
            MARRAY_ASSERT(len_[dim] > 0);

            std::array<len_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, NDim-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, NDim-dim-1, stride.begin()+dim);

            return {len, data_, stride};
        }

        template <int Dim=0, int N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cback() const
        {
            return const_cast<marray_base&>(*this).back<Dim>();
        }

        template <int Dim=0, int N=NDim>
        detail::enable_if_t<N==1, cref>
        back() const
        {
            return const_cast<marray_base&>(*this).back<Dim>();
        }

        template <int Dim=0, int N=NDim>
        detail::enable_if_t<N==1, reference>
        back()
        {
            static_assert(Dim == 0, "Dim out of range");
            MARRAY_ASSERT(len_[0] > 0);
            return data_[(len_[0]-1)*stride_[0]];
        }

        template <int N=NDim>
        detail::enable_if_t<N==1, const_reference>
        cback(int dim) const
        {
            return const_cast<marray_base&>(*this).back(dim);
        }

        template <int N=NDim>
        detail::enable_if_t<N==1, cref>
        back(int dim) const
        {
            return const_cast<marray_base&>(*this).back(dim);
        }

        template <int N=NDim>
        detail::enable_if_t<N==1, reference>
        back(int dim)
        {
            (void)dim;
            MARRAY_ASSERT(dim == 0);
            return back();
        }

        template <int Dim, int N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cback() const
        {
            return const_cast<marray_base&>(*this).back<Dim>();
        }

        template <int Dim, int N=NDim>
        detail::enable_if_t<N!=1, marray_view<ctype, NDim-1>>
        back() const
        {
            return const_cast<marray_base&>(*this).back<Dim>();
        }

        template <int Dim, int N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        back()
        {
            return back<N>(Dim);
        }

        template <int N=NDim>
        detail::enable_if_t<N!=1, marray_view<const Type, NDim-1>>
        cback(int dim) const
        {
            return const_cast<marray_base&>(*this).back(dim);
        }

        template <int N=NDim>
        detail::enable_if_t<N!=1, marray_view<ctype, NDim-1>>
        back(int dim) const
        {
            return const_cast<marray_base&>(*this).back(dim);
        }

        template <int N=NDim>
        detail::enable_if_t<N!=1, marray_view<Type, NDim-1>>
        back(int dim)
        {
            MARRAY_ASSERT(dim >= 0 && dim < NDim);
            MARRAY_ASSERT(len_[dim] > 0);

            std::array<len_type, NDim-1> len;
            std::array<stride_type, NDim-1> stride;

            std::copy_n(len_.begin(), dim, len.begin());
            std::copy_n(len_.begin()+dim+1, NDim-dim-1, len.begin()+dim);
            std::copy_n(stride_.begin(), dim, stride.begin());
            std::copy_n(stride_.begin()+dim+1, NDim-dim-1, stride.begin()+dim);

            return {len, data_ + (len_[dim]-1)*stride_[dim], stride};
        }

        /***********************************************************************
         *
         * Indexing
         *
         **********************************************************************/

        template <int N=NDim>
        detail::enable_if_t<N==1, cref>
        operator[](len_type i) const
        {
            return const_cast<marray_base&>(*this)[i];
        }

        template <int N=NDim>
        detail::enable_if_t<N==1, reference>
        operator[](len_type i)
        {
            MARRAY_ASSERT(i < len_[0]);
            return data_[i*stride_[0]];
        }

        template <int N=NDim>
        detail::enable_if_t<N!=1, marray_slice<ctype, NDim, 1>>
        operator[](len_type i) const
        {
            MARRAY_ASSERT(i < len_[0]);
            return {*this, i};
        }

        template <int N=NDim>
        detail::enable_if_t<N!=1, marray_slice<Type, NDim, 1>>
        operator[](len_type i)
        {
            MARRAY_ASSERT(i < len_[0]);
            return {*this, i};
        }

        template <typename I>
        marray_slice<ctype, NDim, 1, slice_dim>
        operator[](const range_t<I>& x) const
        {
            MARRAY_ASSERT(x.size() >= 0);
            MARRAY_ASSERT(x.front() >= 0);
            MARRAY_ASSERT(x.front()+x.size() <= len_[0]);
            return {*this, x};
        }

        template <typename I>
        marray_slice<Type, NDim, 1, slice_dim>
        operator[](const range_t<I>& x)
        {
            MARRAY_ASSERT(x.size() >= 0);
            MARRAY_ASSERT(x.front() >= 0);
            MARRAY_ASSERT(x.front()+x.size() <= len_[0]);
            return {*this, x};
        }

        marray_slice<ctype, NDim, 1, slice_dim>
        operator[](all_t) const
        {
            return {*this, range(len_[0])};
        }

        marray_slice<Type, NDim, 1, slice_dim>
        operator[](all_t)
        {
            return {*this, range(len_[0])};
        }

        marray_slice<ctype, NDim, 0, bcast_dim>
        operator[](bcast_t) const
        {
            return {*this, slice::bcast};
        }

        marray_slice<Type, NDim, 0, bcast_dim>
        operator[](bcast_t)
        {
            return {*this, slice::bcast, len_[0]};
        }

        cref operator()(detail::array_1d<len_type> idx_) const
        {
            return const_cast<marray_base&>(*this)(idx_);
        }

        reference operator()(detail::array_1d<len_type> idx_)
        {
            MARRAY_ASSERT(idx_.size() == dimension());

            std::array<len_type,NDim> idx;
            idx_.slurp(idx);

            auto ptr = data();
            for (auto i : range(NDim))
                ptr += idx[i]*stride(i);

            return *ptr;
        }

        template <typename Arg, typename=
            detail::enable_if_t<detail::is_index_or_slice<Arg>::value>>
        auto operator()(Arg&& arg) const ->
        decltype((*this)[std::forward<Arg>(arg)])
        {
            return (*this)[std::forward<Arg>(arg)];
        }

        template <typename Arg, typename=
            detail::enable_if_t<detail::is_index_or_slice<Arg>::value>>
        auto operator()(Arg&& arg) ->
        decltype((*this)[std::forward<Arg>(arg)])
        {
            return (*this)[std::forward<Arg>(arg)];
        }

        template <typename Arg, typename... Args, typename=
            detail::enable_if_t<sizeof...(Args) &&
                detail::are_indices_or_slices<Arg, Args...>::value>>
        auto operator()(Arg&& arg, Args&&... args) const ->
        decltype((*this)[std::forward<Arg>(arg)](std::forward<Args>(args)...))
        {
            return (*this)[std::forward<Arg>(arg)](std::forward<Args>(args)...);
        }

        template <typename Arg, typename... Args, typename=
            detail::enable_if_t<sizeof...(Args) &&
                detail::are_indices_or_slices<Arg, Args...>::value>>
        auto operator()(Arg&& arg, Args&&... args) ->
        decltype((*this)[std::forward<Arg>(arg)](std::forward<Args>(args)...))
        {
            return (*this)[std::forward<Arg>(arg)](std::forward<Args>(args)...);
        }

        /***********************************************************************
         *
         * Iteration
         *
         **********************************************************************/

        template <typename Func>
        void for_each_element(Func&& f) const
        {
            for_each_element<cptr>(std::forward<Func>(f),
                                   detail::static_range<int, NDim>{});
        }

        template <typename Func>
        void for_each_element(Func&& f)
        {
            for_each_element<pointer>(std::forward<Func>(f),
                                      detail::static_range<int, NDim>{});
        }

        /***********************************************************************
         *
         * Basic getters
         *
         **********************************************************************/

        const_pointer cdata() const
        {
            return const_cast<marray_base&>(*this).data();
        }

        cptr data() const
        {
            return const_cast<marray_base&>(*this).data();
        }

        pointer data()
        {
            return data_;
        }

        template <typename=void, int N=NDim, typename=detail::enable_if_t<N==1>>
        len_type length() const
        {
            return len_[0];
        }

        template <int Dim>
        len_type length() const
        {
            static_assert(Dim >= 0 && Dim < NDim, "Dim out of range");
            return len_[Dim];
        }

        len_type length(int dim) const
        {
            MARRAY_ASSERT(dim >= 0 && dim < NDim);
            return len_[dim];
        }

        const std::array<len_type, NDim>& lengths() const
        {
            return len_;
        }

        template <typename=void, int N=NDim, typename=detail::enable_if_t<N==1>>
        stride_type stride() const
        {
            return stride_[0];
        }

        template <int Dim>
        stride_type stride() const
        {
            static_assert(Dim >= 0 && Dim < NDim, "Dim out of range");
            return stride_[Dim];
        }

        stride_type stride(int dim) const
        {
            MARRAY_ASSERT(dim >= 0 && dim < NDim);
            return stride_[dim];
        }

        const std::array<stride_type, NDim>& strides() const
        {
            return stride_;
        }

        static constexpr auto dimension()
        {
            return NDim;
        }

        friend std::ostream& operator<<(std::ostream& os, const marray_base& x)
        {
            for (auto i : range(NDim-1))
                os << std::string(i, ' ') << "{\n";

            std::array<len_type,NDim-1> idx = {};
            auto data = x.data_;

            for (bool done = false;!done;)
            {
                os << std::string(NDim-1, ' ') << '{';
                auto n = x.len_[NDim-1];
                if (n > 0)
                {
                    for (auto i : range(n-1))
                        os << data[i*x.stride_[NDim-1]] << ", ";
                    os << data[(n-1)*x.stride_[NDim-1]];
                }
                os << "}";

                for (auto i : reversed_range(NDim-1))
                {
                    idx[i]++;
                    data += x.stride_[i];

                    if (idx[i] >= x.len_[i])
                    {
                        data -= idx[i]*x.stride_[i];
                        idx[i] = 0;
                        os << "\n" << std::string(i, ' ') << '}';
                        if (i == 0) done = true;
                    }
                    else
                    {
                        os << ",\n";
                        for (auto j : range(i+1,NDim-1))
                            os << std::string(j, ' ') << "{\n";
                        break;
                    }
                }

                if (NDim == 1) break;
            }

            return os;
        }
};

}

#endif
