#ifndef MICROHHC_CUDA_TILING_H
#define MICROHHC_CUDA_TILING_H

#ifndef __CUDACC_RTC__
#include "grid.h"
#endif

#ifdef __CUDACC__
#define CUDA_DEVICE __device__ __forceinline__
#define CUDA_HOST_DEVICE __host__ CUDA_DEVICE
#define CUDA_KERNEL __global__
#define CUDA_KERNEL_BOUNDS(threads_per_block, blocks_per_sm) CUDA_KERNEL __launch_bounds__(threads_per_block, blocks_per_sm)
#define CUDA_ASSUME(expr) __builtin_assume(expr)
#define CUDA_EXPECT(expr) __builtin_expect(expr)
#else
#define CUDA_HOST_DEVICE
#define CUDA_ASSUME(expr) do {} while ()
#define CUDA_EXPECT(expr) do {} while ()
#endif

struct DynBlockSize {
    CUDA_DEVICE
    static constexpr int get(int axis) {
        if (axis == 0) return blockDim.x;
        if (axis == 1) return blockDim.y;
        if (axis == 2) return blockDim.z;
        return 1;
    }
};

template <unsigned int X, unsigned int Y, unsigned int Z>
struct StaticBlockSize {
    static_assert(X > 0 && Y > 0 && Z > 0, "invalid block size");

    CUDA_DEVICE
    static constexpr int get(int axis) {
        if (axis == 0) return X;
        if (axis == 1) return Y;
        if (axis == 2) return Z;
        return 1;
    }
};

template <typename BlockSize = DynBlockSize>
struct TilingStrategy {
    CUDA_HOST_DEVICE
    constexpr TilingStrategy(
            unsigned int tile_factor_x,
            unsigned int tile_factor_y,
            unsigned int tile_factor_z,
            unsigned int unroll_factor_x,
            unsigned int unroll_factor_y,
            unsigned int unroll_factor_z,
            bool tile_contiguous_x,
            bool tile_contiguous_y,
            bool tile_contiguous_z,
            unsigned int blocks_per_sm_
    ):
            tile_factor_ {tile_factor_x, tile_factor_y, tile_factor_z},
            unroll_factor_ {unroll_factor_x, unroll_factor_y, unroll_factor_z},
            tile_contiguous_ {tile_contiguous_x, tile_contiguous_y, tile_contiguous_z},
            blocks_per_sm_ {blocks_per_sm_} {}

    CUDA_HOST_DEVICE
    constexpr TilingStrategy(): TilingStrategy(
            32, 4, 1,
            1, 1, 1,
            1, 1, 1,
            false, false, false,
            1) {}

    CUDA_HOST_DEVICE
    constexpr unsigned int block_size(size_t axis) const {
        if (axis >= 3) return 1;
        return BlockSize::get(axis);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_factor(size_t axis) const {
        if (axis >= 3 || tile_factor_[axis] == 0) return 1;
        return tile_factor_[axis];
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_size(size_t axis) const {
        return tile_factor(axis) * block_size(axis);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_contiguous(size_t axis) const {
        if (axis >= 3) return false;
        return tile_contiguous_[axis];
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_contiguous_x() const {
        return tile_contiguous(0);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_contiguous_y() const {
        return tile_contiguous(1);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_contiguous_z() const {
        return tile_contiguous(2);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int unroll_factor(size_t axis) const {
        if (axis >= 3) return 1;
        return unroll_factor_[axis] > 0 ? unroll_factor_[axis] : tile_factor(axis);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int threads_per_block() const {
        return block_size(0) * block_size(1) * block_size(2);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int items_per_thread() const {
        return tile_factor(0) * tile_factor(1) * tile_factor(2);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int items_per_block() const {
        return threads_per_block() * items_per_thread();
    }

    CUDA_HOST_DEVICE
    constexpr dim3 block_size() const {
        return dim3(block_size(0), block_size(1), block_size(2));
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int block_size_x() const {
        return block_size(0);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int block_size_y() const {
        return block_size(1);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int block_size_z() const {
        return block_size(2);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_factor_x() const {
        return tile_factor(0);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_factor_y() const {
        return tile_factor(1);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_factor_z() const {
        return tile_factor(2);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_size_x() const {
        return tile_size(0);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_size_y() const {
        return tile_size(1);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int tile_size_z() const {
        return tile_size(2);
    }

    CUDA_HOST_DEVICE
    constexpr unsigned int blocks_per_sm() const {
        return blocks_per_sm_;
    }

private:
    unsigned int tile_factor_[3];
    unsigned int unroll_factor_[3];
    bool tile_contiguous_[3];
    unsigned int blocks_per_sm_;
};

template <
        unsigned int block_size_x_,
        unsigned int block_size_y_,
        unsigned int block_size_z_,
        unsigned int tile_factor_x_,
        unsigned int tile_factor_y_,
        unsigned int tile_factor_z_,
        unsigned int unroll_factor_x,
        unsigned int unroll_factor_y,
        unsigned int unroll_factor_z,
        bool tile_contiguous_x_,
        bool tile_contiguous_y_,
        bool tile_contiguous_z_,
        unsigned int blocks_sm_ = 1
>
struct StaticTilingStrategy: TilingStrategy<StaticBlockSize<
        block_size_x_,
        block_size_y_,
        block_size_z_
    >> {
    using base_type = TilingStrategy<StaticBlockSize<block_size_x_, block_size_y_, block_size_z_>>;


    CUDA_HOST_DEVICE
    constexpr StaticTilingStrategy(): base_type(
            tile_factor_x_,
            tile_factor_y_,
            tile_factor_z_,
            unroll_factor_x,
            unroll_factor_y,
            unroll_factor_z,
            tile_contiguous_x_,
            tile_contiguous_y_,
            tile_contiguous_z_,
            blocks_sm_
    ) {}
};

using DefaultTilingStrategy = StaticTilingStrategy<128, 1, 1, 1, 1, 1, 1, 1, 1, false, false, false>;

struct Grid_layout
{
    const int istart;
    const int jstart;
    const int kstart;
    const int iend;
    const int jend;
    const int kend;
    const int ii;
    const int jj;
    const int kk;

#ifndef __CUDACC_RTC__
    template <typename TF>
    static Grid_layout from_grid_data(const Grid_data<TF>& gd)
    {
        return {
            .istart = gd.istart,
            .jstart = gd.jstart,
            .kstart = gd.kstart,
            .iend = gd.iend,
            .jend = gd.jend,
            .kend = gd.kend,
            .ii = 1,
            .jj = gd.jcells,
            .kk = gd.ijcells
        };
    }
#endif

    CUDA_HOST_DEVICE
    int operator()(int i, int j, int k) const
    {
        return ii * i + jj * j + kk * k;
    }
};

namespace levels
{
    struct Interior
    {
        constexpr static int MAX_VALUE = ~(1 << 31);

        CUDA_HOST_DEVICE
        int distance_to_start() const
        {
            return MAX_VALUE;
        }

        CUDA_HOST_DEVICE
        int distance_to_end() const
        {
            return MAX_VALUE;
        }
    };

    struct General
    {
        CUDA_HOST_DEVICE
        General(int k, int kstart, int kend)
        {
            dist_start_ = k - kstart;
            dist_end_ = kend - k - 1;
        }


        CUDA_HOST_DEVICE
        int distance_to_start() const
        {
            return dist_start_;
        }

        CUDA_HOST_DEVICE
        int distance_to_end() const
        {
            return dist_end_;
        }

    private:
        int dist_start_;
        int dist_end_;
    };
}

template <typename F, typename Tiling = DefaultTilingStrategy, size_t BorderSize = 3>
struct TilingKernel
{
    static constexpr Tiling tiling {};
    static constexpr size_t border_size = BorderSize;
    static constexpr dim3 block_size = tiling.block_size();
    static constexpr int blocks_per_sm = tiling.blocks_per_sm();

    static_assert(tiling.tile_size_x() > 0, "invalid tile_size_x");
    static_assert(tiling.tile_size_y() > 0, "invalid tile_size_y");
    static_assert(tiling.tile_size_z() > 0, "invalid tile_size_z");

    CUDA_HOST_DEVICE
    TilingKernel(F fun = {}) : fun_(fun) {}

    template <typename... Args>
    CUDA_DEVICE
    void process_cta_layer(const Grid_layout& grid, F fun, Args... args)
    {
        const int thread_idx_x = tiling.block_size_x() > 1 ? threadIdx.x : 0;
        const int xoffset = grid.istart + blockIdx.x * tiling.tile_size_x() +
                            (tiling.tile_contiguous_x() ? thread_idx_x * tiling.tile_factor_x() : thread_idx_x);

        if (tiling.block_size_x() > 1 && tiling.tile_factor_x() == 1 && xoffset >= grid.iend) return;

        const int thread_idx_y = tiling.block_size_y() > 1 ? threadIdx.y : 0;
        const int yoffset = grid.jstart + blockIdx.y * tiling.tile_size_y() +
                            (tiling.tile_contiguous_y() ? thread_idx_y * tiling.tile_factor_y() : thread_idx_y);

        if (tiling.block_size_y() > 1 && tiling.tile_factor_y() == 1 && yoffset >= grid.jend) return;

#pragma unroll(tiling.unroll_factor(1))
        for (int dy = 0; dy < tiling.tile_factor_y(); dy++)
        {
            const int j = yoffset + dy * (tiling.tile_contiguous_y() ? 1 : tiling.block_size_y());
            if (tiling.tile_factor_y() > 1 && j >= grid.jend) break;

#pragma unroll(tiling.unroll_factor(0))
            for (int dx = 0; dx < tiling.tile_factor_x(); dx++)
            {
                const int i = xoffset + dx * (tiling.tile_contiguous_x() ? 1 : tiling.block_size_x());
                if (tiling.tile_factor_x() > 1 && i >= grid.iend) break;

                fun_(grid, i, j, args...);
            }
        }
    }

    template <typename... Args>
    CUDA_DEVICE
    void process_cta(const Grid_layout& grid, F fun, Args... args)
    {
        const int thread_idx_z = tiling.block_size_z() > 1 ? threadIdx.z : 0;
        const int block_idx_z = blockIdx.z;
        const int zoffset = grid.kstart + block_idx_z * tiling.tile_size_z() +
                thread_idx_z * (tiling.tile_contiguous_z() ?  tiling.tile_factor_z() : 1);
        const int zstep = (tiling.tile_contiguous_z() ? 1 : tiling.block_size_z());

        const int zlow = zoffset;
        const int zhigh = zoffset + (tiling.tile_factor_z() - 1) * zstep;

        if (border_size != 0 && zlow >= grid.kstart + border_size && zhigh <= grid.kend - border_size + 1)
        {
#pragma unroll (tiling.unroll_factor(2))
            for (int dz = 0; dz < tiling.tile_factor_z(); dz++)
            {
                const int k = zoffset + dz * zstep;
                levels::Interior level;
                process_cta_layer(grid, fun, k, level, args...);
            }
        }
        else
        {
#pragma unroll (tiling.unroll_factor(2))
            for (int dz = 0; dz < tiling.tile_factor_z(); dz++)
            {
                const int k = zoffset + dz * zstep;
                if (tiling.tile_size_z() > 1 && k >= grid.kend) break;

                levels::General level {k, grid.kstart, grid.kend};
                process_cta_layer(grid, fun, k, level, args...);
            }
        }
    }

    template <typename... Args>
    CUDA_DEVICE
    void operator()(const Grid_layout& grid, Args... args)
    {
        process_cta(grid, F{}, args...);
    }

private:
    F fun_;
};


template <typename F, typename... Args>
//CUDA_KERNEL_BOUNDS(F::block_size.x * F::block_size.y * F::block_size.z, F::blocks_per_sm)
__global__
void kernel_wrapper(Args... args)
{
    F{}(args...);
}

/// TODO: move to better location.
struct SourceLocation
{

    CUDA_HOST_DEVICE
    constexpr SourceLocation(const char* file_name, int line) : file_name_(file_name), line_(line) {}

    CUDA_HOST_DEVICE
    const char* file_name() const
    {
        return file_name_;
    }

    CUDA_HOST_DEVICE
    int lineno() const
    {
        return line_;
    }

private:
    const char* const file_name_;
    const int line_;
};

#define DEFINE_KERNEL_CONSTANTS \
    static constexpr SourceLocation source_location = SourceLocation(__FILE__, __LINE__);


#endif //MICROHHC_CUDA_TILING_H
