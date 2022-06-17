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


/// Represents arbitrary memory allocated on CUDA device.
struct cuda_raw_buffer
{
    explicit cuda_raw_buffer(size_t size_in_bytes = 0);
    cuda_raw_buffer(cuda_raw_buffer&& that) noexcept;
    void resize(size_t size_in_bytes);

    cuda_raw_buffer& operator=(cuda_raw_buffer&&);

    cuda_raw_buffer(const cuda_raw_buffer&) = delete;
    cuda_raw_buffer& operator=(const cuda_raw_buffer&) = delete;

    void *get()
    {
        return ptr_;
    }

    const void* get() const
    {
        return ptr_;
    }

private:
    void *ptr_ = nullptr;
};

/// Templated wrapper around cuda_raw_buffer that adds type information.
template <typename T>
struct cuda_buffer
{
    cuda_buffer(size_t size = 0)
    {
        resize(size);
    }

    cuda_buffer(cuda_buffer&& that)
    {
        std::swap(that.buffer_, buffer_);
        std::swap(that.size_, size_);
    }

    cuda_buffer(const cuda_buffer&) = delete;

    ~cuda_buffer()
    {
        free();
    }

    T* get()
    {
        return static_cast<T*>(buffer_.get());
    }

    const T* get() const
    {
        return static_cast<const T*>(buffer_.get());
    }

    // FIXME: should these operator overloads be explicit?
    operator T*()
    {
        return get();
    }

    operator const T*() const
    {
        return get();
    }

    void resize(size_t new_size)
    {
        size_ = 0;  // Set size to zero since buffer_.resize() might throw
        buffer_.resize(new_size * sizeof(T));
        size_ = new_size;
    }

    void free()
    {
        resize(0);
    }

    size_t size() const
    {
        return size_;
    }

    size_t size_in_bytes() const
    {
        return size() * sizeof(T);
    }

private:
    size_t size_ = 0;
    cuda_raw_buffer buffer_;
};

#endif //CUDA_BUFFER_H
