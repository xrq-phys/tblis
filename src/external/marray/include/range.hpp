#ifndef _MARRAY_RANGE_HPP_
#define _MARRAY_RANGE_HPP_

#include "utility.hpp"

namespace MArray
{

namespace detail
{

template <typename T, typename=void>
struct underlying_type_if
{
    typedef T type;
};

template <typename T>
struct underlying_type_if<T, detail::enable_if_t<std::is_enum<T>::value>>
{
    typedef typename std::underlying_type<T>::type type;
};

template <typename... Ts> struct are_numeric;

template <> struct are_numeric<>
: std::integral_constant<bool,true> {};

template <typename T, typename... Ts> struct are_numeric<T, Ts...>
: std::integral_constant<bool, (std::is_integral<T>::value ||
                                std::is_enum<T>::value) &&
                               are_numeric<Ts...>::value> {};

template <typename... Ts> using enable_if_numeric =
    std::enable_if_t<are_numeric<Ts...>::value>;

}

template <typename T>
class range_t
{
    static_assert(std::is_integral<T>::value, "The type must be integral.");

    protected:
        T from_;
        T to_;
        T delta_;

        typedef T value_type;
        typedef T size_type;

    public:
        class iterator : std::iterator<std::random_access_iterator_tag,T>
        {
            protected:
                T val_;
                T delta_;

            public:
                using typename std::iterator<std::random_access_iterator_tag,T>::iterator_category;
                using typename std::iterator<std::random_access_iterator_tag,T>::value_type;
                using typename std::iterator<std::random_access_iterator_tag,T>::difference_type;
                using typename std::iterator<std::random_access_iterator_tag,T>::pointer;
                using typename std::iterator<std::random_access_iterator_tag,T>::reference;

                constexpr iterator() : val_(0), delta_(0) {}

                constexpr iterator(T val, T delta) : val_(val), delta_(delta) {}

                bool operator==(const iterator& other) const
                {
                    return val_ == other.val_ && delta_ == other.delta_;
                }

                bool operator!=(const iterator& other) const
                {
                    return val_ != other.val_ || delta_ != other.delta_;
                }

                value_type operator*() const
                {
                    return val_;
                }

                iterator& operator++()
                {
                    val_ += delta_;
                    return *this;
                }

                iterator operator++(int)
                {
                    iterator old(*this);
                    val_ += delta_;
                    return old;
                }

                iterator& operator--()
                {
                    val_ -= delta_;
                    return *this;
                }

                iterator operator--(int)
                {
                    iterator old(*this);
                    val_ -= delta_;
                    return old;
                }

                iterator& operator+=(difference_type n)
                {
                    val_ += n*delta_;
                    return *this;
                }

                iterator operator+(difference_type n) const
                {
                    return iterator(val_+n*delta_, delta_);
                }

                friend iterator operator+(difference_type n, const iterator& i)
                {
                    return iterator(i.val_+n*i.delta_, i.delta_);
                }

                iterator& operator-=(difference_type n)
                {
                    val_ -= n*delta_;
                    return *this;
                }

                iterator operator-(difference_type n) const
                {
                    return iterator(val_-n*delta_, delta_);
                }

                difference_type operator-(const iterator& other) const
                {
                    return (val_-other.val_)/delta_;
                }

                bool operator<(const iterator& other) const
                {
                    return val_ < other.val_;
                }

                bool operator<=(const iterator& other) const
                {
                    return val_ <= other.val_;
                }

                bool operator>(const iterator& other) const
                {
                    return val_ > other.val_;
                }

                bool operator>=(const iterator& other) const
                {
                    return val_ >= other.val_;
                }

                value_type operator[](difference_type n) const
                {
                    return val_+n*delta_;
                }

                friend void swap(iterator& a, iterator& b)
                {
                    using std::swap;
                    swap(a.val_, b.val_);
                    swap(a.delta_, b.delta_);
                }
        };

        constexpr range_t()
        : from_(0), to_(0), delta_(0) {}

        constexpr range_t(T to)
        : from_(0), to_(to), delta_(1) {}

        constexpr range_t(T from, T to)
        : from_(from), to_(to), delta_(1) {}

        constexpr range_t(T from, T to, T delta)
        : from_(from), delta_(delta)
        {
            if (delta > 0)
                to_ = from+((to-from+delta-1)/delta)*delta;
            else if (delta < 0)
                to_ = from+((to-from+delta+1)/delta)*delta;
        }

        range_t(const range_t&) = default;

        range_t(range_t&&) = default;

        range_t& operator=(const range_t&) = default;

        range_t& operator=(range_t&&) = default;

        value_type step() const
        {
            return delta_;
        }

        size_type size() const
        {
            return (delta_ == 0 ? std::numeric_limits<size_type>::max() :
                    (to_-from_)/delta_);
        }

        iterator begin() const
        {
            return iterator(from_, delta_);
        }

        iterator end() const
        {
            return iterator(to_, delta_);
        }

        value_type front() const
        {
            return from_;
        }

        value_type back() const
        {
            return to_-delta_;
        }

        value_type operator[](size_type n) const
        {
            return from_+n*delta_;
        }

        template <typename U, typename=
            typename std::enable_if<detail::is_container<U>::value>::type>
        operator U() const
        {
            return U(begin(), end());
        }

        range_t& operator+=(T shift)
        {
            from_ += shift;
            to_ += shift;
            return *this;
        }

        range_t& operator-=(T shift)
        {
            from_ -= shift;
            to_ -= shift;
            return *this;
        }

        range_t operator+(T shift)
        {
            range_t shifted(*this);
            shifted += shift;
            return shifted;
        }

        range_t operator-(T shift)
        {
            range_t shifted(*this);
            shifted -= shift;
            return shifted;
        }

        friend range_t operator+(T shift, const range_t& other)
        {
            return other + shift;
        }

        friend range_t operator-(T shift, const range_t& other)
        {
            range_t shifted(other);
            shifted.from_ = shift - shifted.from_;
            shifted.to_ -= shift - shifted.to_;
            shifted.delta_ = -shifted.delta_;
            return shifted;
        }
};

template <typename T, typename=
    detail::enable_if_numeric<T>>
auto range(T to)
{
    typedef typename detail::underlying_type_if<T>::type U;
    return range_t<U>{U(to)};
}

template <typename T, typename U, typename=
    detail::enable_if_numeric<T,U>>
auto rangeN(T from, U N)
{
    typedef decltype(std::declval<T>() + std::declval<U>()) V0;
    typedef typename detail::underlying_type_if<V0>::type V;
    return range_t<V>{V(from), V(from+N)};
}

template <typename T, typename U, typename=
    detail::enable_if_numeric<T,U>>
auto range(T from, U to)
{
    typedef decltype(std::declval<T>() + std::declval<U>()) V0;
    typedef typename detail::underlying_type_if<V0>::type V;
    if ((V)to < (V)from)
        to = from;
    return range_t<V>{(V)from, (V)to};
}

template <typename T, typename U, typename V, typename=
    detail::enable_if_numeric<T,U,V>>
auto range(T from, U to, V delta)
{
    typedef decltype(std::declval<T>() + std::declval<U>() + std::declval<V>()) W0;
    typedef typename detail::underlying_type_if<W0>::type W;
    if ((W() < (W)delta && (W)to < (W)from) ||
        ((W)delta < W() && (W)from < (W)to))
        to = from;
    return range_t<W>{(W)from, (W)to, (W)delta};
}

template <typename T, typename=
    detail::enable_if_numeric<T>>
auto reversed_range(T to)
{
    return range(to-1, -1, -1);
}

template <typename T, typename U, typename=
    detail::enable_if_numeric<T,U>>
auto reversed_rangeN(T from, U N)
{
    return range(from-1, from-N-1, -1);
}

template <typename T, typename U, typename=
    detail::enable_if_numeric<T,U>>
auto reversed_range(T from, U to)
{
    return range(to-1, from-1, -1);
}

}

#endif
