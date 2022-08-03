/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef CUDA_BUFFER_H
#define CUDA_BUFFER_H

#include <stdexcept>
#include <vector>
#include <limits>

template <typename Derived, typename T, bool Owned>
struct cuda_buffer_base;

template <typename T>
struct cuda_span;

template <typename T>
struct cuda_vector;

/**
 * Copy bytes between two memory location using `cudaMemcpy`. The copy is performed without any type checking.
 *
 * @param src The source location.
 * @param dst The destination location.
 * @param nbytes The number of bytes.
 */
void cuda_raw_copy(const void* src, void* dst, size_t nbytes);

/**
 * Copy `num_elements` elements from `src` to `dst` using `cudaMemcpy`.
 */
template <typename T, typename = typename std::enable_if<std::is_trivially_copyable<T>::value>::type>
void cuda_copy(const T* src, T* dst, size_t num_elements)
{
    cuda_raw_copy(
            static_cast<const void*>(src),
            static_cast<void*>(dst),
            num_elements * sizeof(T));
}

template <typename T> struct cuda_copy_location;
template <typename T> struct cuda_copy_location<T&>: cuda_copy_location<T> {};
template <typename T> struct cuda_copy_location<T&&>: cuda_copy_location<T> {};

template <typename T> struct cuda_copy_location<cuda_span<T>> { using type = T; };
template <typename T> struct cuda_copy_location<const cuda_span<T>> { using type = T; };

template <typename T> struct cuda_copy_location<cuda_vector<T>> { using type = T; };
template <typename T> struct cuda_copy_location<const cuda_vector<T>> { using type = const T; };

template <typename T> struct cuda_copy_location<std::vector<T>> { using type = T; };
template <typename T> struct cuda_copy_location<const std::vector<T>> { using type = const T; };

template <typename A, typename B>
struct is_cuda_copy_compatible: std::integral_constant<bool,
        std::is_trivially_copyable<A>::value &&
        std::is_trivially_copyable<B>::value &&
        std::is_same<typename std::remove_cv<A>::type, B>::value
        > {};

/**
 * Copy elements from `src` to `dst` using `cudaMemcpy`. An exception is thrown if the sizes of both buffers do not
 * match.
 */
template <typename Src, typename Dst, typename = typename std::enable_if<is_cuda_copy_compatible<
    typename cuda_copy_location<Src>::type,
    typename cuda_copy_location<Dst>::type
>::value>::type>
void cuda_copy(Src&& src, Dst&& dst)
{
    using src_type = typename cuda_copy_location<Src>::type;
    using dst_type = typename cuda_copy_location<Dst>::type;

    cuda_span<const src_type> src_view = {src.data(), src.size()};
    cuda_span<dst_type> dst_view = {dst.data(), dst.size()};

    if (src_view.size_in_bytes() != dst_view.size_in_bytes())
    {
        throw std::runtime_error("buffer sizes do not match");
    }

    cuda_raw_copy(
            static_cast<const void*>(src_view.data()),
            static_cast<void*>(dst_view.data()),
            src.size_in_bytes());
}

/**
 * Represents memory buffer allocated on a CUDA device using `cudaMalloc`. Memory is automatically freed on
 * destruction.
 */
struct cuda_raw_buffer
{
    explicit cuda_raw_buffer(size_t size_in_bytes = 0);
    cuda_raw_buffer(cuda_raw_buffer&& that) noexcept;
    cuda_raw_buffer& operator=(cuda_raw_buffer&&) noexcept;
    cuda_raw_buffer(const cuda_raw_buffer&) = delete;
    cuda_raw_buffer& operator=(const cuda_raw_buffer&) = delete;

    /**
     * Reallocate memory. This essentially calls `cudaFree` on the old pointer, followed by `cudaMalloc` on the new
     * pointer. If new allocation fails, this method throws an exception and the allocation is set to `nullptr`.
     *
     * @param size_in_bytes
     */
    void reallocate(size_t size_in_bytes);

    /**
     * Pointer to the data
     */
    void *data() noexcept
    {
        return ptr_;
    }

    /**
     * Pointer to the data.
     */
    const void* data() const noexcept
    {
        return ptr_;
    }

private:
    void *ptr_ = nullptr;
};


/**
 * Represents a non-owned array of elements of type `T` in GPU memory. This is essentially just a wrapper around a
 * (ptr, size) pair.
 */
template <typename T>
struct cuda_span: cuda_buffer_base<cuda_span<T>, T, false>
{
    cuda_span(T* ptr, size_t size): ptr_(ptr), size_(size) {}
    cuda_span() = default;
    cuda_span(const cuda_span&) = default;
    cuda_span(cuda_span&&) noexcept = default;

    T* data() const noexcept
    {
        return ptr_;
    }

    size_t size() const noexcept
    {
        return size_;
    }

    operator T*() const noexcept
    {
        return data();
    }

    operator const T*() const noexcept
    {
        return data();
    }

    operator cuda_span<const T>() const noexcept
    {
        return {data(), size()};
    }

private:
    T* ptr_ = nullptr;
    size_t size_ = 0;
};

/**
 * Represents a non-owned array of elements of type `T` in GPU memory. This is essentially just a wrapper around a
 * (ptr, size) pair. Specialization for `const T~.
 */
template <typename T>
struct cuda_span<const T>: cuda_buffer_base<cuda_span<const T>, const T, false>
{
    cuda_span() = default;
    cuda_span(const T* ptr, size_t size): ptr_(ptr), size_(size) {}
    cuda_span(const cuda_span&) = default;
    cuda_span(cuda_span&&) noexcept = default;

    const T* data() const noexcept
    {
        return ptr_;
    }

    size_t size() const noexcept
    {
        return size_;
    }

    operator const T*() const noexcept
    {
        return data();
    }

private:
    const T* ptr_ = nullptr;
    size_t size_ = 0;
};


/**
 * Represents RAII wrapper around an array of `T` allocated on CUDA device using `cudaMalloc`.
 */
template <typename T>
struct cuda_vector: cuda_buffer_base<cuda_vector<T>, T, true>
{
    cuda_vector(size_t size = 0)
    {
        allocate(size);
    }

    cuda_vector(const std::vector<T>& data)
    {
        from_vector(data);
    }

    cuda_vector(std::initializer_list<T> data)
    {
        from_vector(data);
    }

    cuda_vector(cuda_vector&& that) noexcept
    {
        *this = std::move(that);
    }

    cuda_vector& operator=(cuda_vector&& that) noexcept
    {
        std::swap(that.buffer_, buffer_);
        std::swap(that.size_, size_);
        return *this;
    }

    // No copy constructor/assignment
    cuda_vector& operator=(const cuda_vector&) = delete;
    cuda_vector(const cuda_vector&) = delete;

    ~cuda_vector()
    {
        free();
    }

    /**
     * Pointer to the data
     */
    T* data() noexcept
    {
        return static_cast<T*>(buffer_.data());
    }

    /**
     * Pointer to the data
     */
    const T* data() const noexcept
    {
        return static_cast<const T*>(buffer_.data());
    }

    /**
     * The size of this span given as the number of elements.
     */
    size_t size() const noexcept
    {
        return size_;
    }

    operator T*() noexcept
    {
        return data();
    }

    operator const T*() const noexcept
    {
        return data();
    }

    operator cuda_span<T>() noexcept
    {
        return {data(), size()};
    }

    operator cuda_span<const T>() const noexcept
    {
        return {data(), size()};
    }

    /**
     * Allocate memory for `new_size` elements. Existing data in this vector is not preserved.
     *
     * @param new_size New size of the vector in number of elements.
     */
    void allocate(size_t new_size)
    {
        // resize might throw, so we temporarily set the size to zero
        size_ = 0;
        buffer_.reallocate(new_size * sizeof(T));
        size_ = new_size;
    }

    /**
     * Release buffer of this vector.
     */
    void free()
    {
        allocate(0);
    }

    /**
     * Allocate memory for this vector. Existing data in this vector is copied to the new memory location.
     *
     * @param new_size New size of the vector in number of elements.
     */
    void resize(size_t new_size)
    {
        if (new_size != size_)
        {
            cuda_vector<T> new_buffer(new_size);
            size_t noverlap = std::min(new_size, size_);
            cuda_copy(this->subspan(0, noverlap), new_buffer.subspan(0, noverlap));
            *this = std::move(new_buffer);
        }
    }

    /**
     * Resize this `cuda_vector` and copy data from the given host vector.
     */
    void from_vector(const std::vector<T>& input) const
    {
        allocate(input.size());
        this->copy_from(input);
    }

private:
    size_t size_ = 0;
    cuda_raw_buffer buffer_;
};


template <typename Derived, typename T, bool Owned>
struct cuda_buffer_base
{
    using value_type = T;
    using const_value_type = typename std::conditional<Owned, const T, T>::type;
    using decay_value_type = typename std::remove_const<T>::type;
    using byte_type = char;

private:
    Derived& derived() noexcept
    {
        return *static_cast<Derived*>(this);
    }

    const Derived& derived() const noexcept
    {
        return *static_cast<const Derived*>(this);
    }

public:
    cuda_span<value_type> view() noexcept
    {
        return cuda_span<value_type>(derived().data(), derived().size());
    }

    cuda_span<const_value_type> view() const noexcept
    {
        return cuda_span<const_value_type>(derived().data(), derived().size());
    }

    cuda_span<const value_type> cview() const noexcept
    {
        return cuda_span<value_type>(derived().data(), derived().size());
    }

private:
    void assert_in_range(size_t offset, size_t count) const
    {
        size_t size = derived().size();

        // We check for two conditions:
        // - offset + count overflows
        // - offset + count exceeds size
        bool overflows = count > std::numeric_limits<size_t>::max() - offset;
        bool out_of_bounds = offset + count > size;

        if (overflows || out_of_bounds)
        {
            // FIXME: Need more detailed error message
            throw std::runtime_error("index out of bounds");
        }
    }

public:
    /**
     * Returns a view `count` elements of this memory allocation starting at offset `offset`.
     */
    cuda_span<value_type> subspan(size_t offset, size_t count)
    {
        assert_in_range(offset, count);
        return cuda_span<value_type>(derived().data() + offset, count);
    }

    /**
     * Returns a view `count` elements of this memory allocation starting at offset `offset`.
     */
    cuda_span<const_value_type> subspan(size_t offset, size_t count) const
    {
        assert_in_range(offset, count);
        return cuda_span<const_value_type>(derived().data() + offset, count);
    }

    /**
     * Size of memory allocation number of bytes.
     */
    size_t size_in_bytes() const noexcept
    {
        return derived().size() * sizeof(T);
    }

    /**
     * Returns a view of the raw bytes in this allocation.
     */
    cuda_span<typename std::conditional<std::is_const<value_type>::value, const byte_type, byte_type>::type>
            view_bytes() noexcept
    {
        return {(char*) derived().data(), size_in_bytes()};
    }

    /**
     * Returns a constant view of the raw bytes in this allocation.
     */
    cuda_span<typename std::conditional<std::is_const<const_value_type>::value, const byte_type, byte_type>::type>
            view_bytes() const noexcept
    {
        return {(char*) derived().data(), size_in_bytes()};
    }

    std::vector<decay_value_type> to_vector() const
    {
        std::vector<decay_value_type> result(derived().size());
        cuda_copy(view(), result);
        return result;
    }
};

#endif