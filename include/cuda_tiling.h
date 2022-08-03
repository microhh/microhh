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
            bool tile_contiguous_z
    ):
            tile_factor_ {tile_factor_x, tile_factor_y, tile_factor_z},
            unroll_factor_ {unroll_factor_x, unroll_factor_y, unroll_factor_z},
            tile_contiguous_ {tile_contiguous_x, tile_contiguous_y, tile_contiguous_z} {}

    CUDA_HOST_DEVICE
    constexpr TilingStrategy(): TilingStrategy(
            32, 4, 1,
            1, 1, 1,
            1, 1, 1,
            false, false, false) {}

    CUDA_HOST_DEVICE
    constexpr unsigned int block_size(size_t axis) const {
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

private:
    unsigned int tile_factor_[3];
    unsigned int unroll_factor_[3];
    bool tile_contiguous_[3];
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
        bool tile_contiguous_z_
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
            tile_contiguous_z_
    ) {}
};

struct DefaultTilingStrategy: TilingStrategy<DynBlockSize> {
    using base_type = TilingStrategy<DynBlockSize>;

    CUDA_HOST_DEVICE
    constexpr DefaultTilingStrategy(): base_type(
            1, 1, 1, 1, 1, 1, false, false, false
    ) {}
};

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

    CUDA_HOST_DEVICE
    bool operator==(const Grid_layout& that) const
    {
        return that.istart == istart &&
                that.jstart == jstart &&
                that.kstart == kstart &&
                that.iend == iend &&
                that.jend == jend &&
                that.kend == kend &&
                that.ii == ii &&
                that.jj == jj &&
                that.kk == kk;
    }

    CUDA_HOST_DEVICE
    bool operator!=(const Grid_layout& that) const
    {
        return !(*this == that);
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

template <typename Tiling, typename F, typename... Args>
CUDA_DEVICE
void cta_execute_tiling_layer(
        const Grid_layout& grid,
        F fun,
        Args&&... args)
{
    static constexpr Tiling tiling = {};

    const int thread_idx_x = tiling.block_size(0) > 1 ? threadIdx.x : 0;
    const int xoffset = grid.istart + blockIdx.x * tiling.tile_size(0) +
                        (tiling.tile_contiguous(0) ? thread_idx_x * tiling.tile_factor(0) : thread_idx_x);
    const int xstep = (tiling.tile_contiguous(0) ? 1 : tiling.block_size(0));

    if (tiling.block_size(0) > 1 && tiling.tile_factor(0) == 1 && xoffset >= grid.iend) return;

    const int thread_idx_y = tiling.block_size(1) > 1 ? threadIdx.y : 0;
    const int yoffset = grid.jstart + blockIdx.y * tiling.tile_size(1) +
                        (tiling.tile_contiguous(1) ? thread_idx_y * tiling.tile_factor(1) : thread_idx_y);
    const int ystep = (tiling.tile_contiguous(1) ? 1 : tiling.block_size(1));

    if (tiling.block_size(1) > 1 && tiling.tile_factor(1) == 1 && yoffset >= grid.jend) return;

#pragma unroll(tiling.unroll_factor(1))
    for (int dy = 0; dy < tiling.tile_factor(1); dy++)
    {
        const int j = yoffset + dy * ystep;
        if (tiling.tile_factor(1) > 1 && j >= grid.jend) break;

#pragma unroll(tiling.unroll_factor(0))
        for (int dx = 0; dx < tiling.tile_factor(0); dx++)
        {
            const int i = xoffset + dx * xstep;
            if (tiling.tile_factor(0) > 1 && i >= grid.iend) break;

            fun(grid, i, j, args...);
        }
    }
}

template <int border_size, typename Tiling, typename F, typename... Args>
CUDA_DEVICE
void cta_execute_tiling_border(
    const Grid_layout& grid,
    F fun,
    Args&&... args
) {
    static constexpr Tiling tiling = {};

    const int thread_idx_z = tiling.block_size(2) > 1 ? threadIdx.z : 0;
    const int block_idx_z = blockIdx.z;
    const int zoffset = grid.kstart + block_idx_z * tiling.tile_size(2) +
                        thread_idx_z * (tiling.tile_contiguous(2) ?  tiling.tile_factor(2) : 1);
    const int zstep = (tiling.tile_contiguous(2) ? 1 : tiling.block_size(2));

    const int zlow = zoffset;
    const int zhigh = zoffset + (tiling.tile_factor(2) - 1) * zstep;

    if (border_size != 0 && zlow >= grid.kstart + border_size && zhigh <= grid.kend - border_size + 1)
    {
#pragma unroll (tiling.unroll_factor(2))
        for (int dz = 0; dz < tiling.tile_factor(2); dz++)
        {
            const int k = zoffset + dz * zstep;
            levels::Interior level;
            cta_execute_tiling_layer<Tiling>(grid, fun, k, level, args...);
        }
    }
    else
    {
#pragma unroll (tiling.unroll_factor(2))
        for (int dz = 0; dz < tiling.tile_factor(2); dz++)
        {
            const int k = zoffset + dz * zstep;
            if (tiling.tile_size(2) > 1 && k >= grid.kend) break;

            levels::General level {k, grid.kstart, grid.kend};
            cta_execute_tiling_layer<Tiling>(grid, fun, k, level, args...);
        }
    }
}
template <typename Tiling = DefaultTilingStrategy, typename F, typename... Args>
CUDA_DEVICE
void cta_execute_tiling(
        const Grid_layout& grid,
        F fun,
        Args&&... args
) {
    cta_execute_tiling_border<0, Tiling>(
        grid, fun, args...);
}


template <typename F, typename... Args>
__global__
void tiling_kernel(Grid_layout grid, Args... args) {
    cta_execute_tiling(grid, F{}, args...);
}

#ifndef __CUDACC_RTC__
/// TODO: move to better location.
struct GridFunctorMeta
{
    CUDA_HOST_DEVICE
    constexpr GridFunctorMeta(const char* file, int line, const char* name) :
            file_(file), line_(line), name_(name) {}

    CUDA_HOST_DEVICE
    const char* name() const
    {
        return name_;
    }

    CUDA_HOST_DEVICE
    const char* file() const
    {
        return file_;
    }

    CUDA_HOST_DEVICE
    int lineno() const
    {
        return line_;
    }

private:
    const char* const file_;
    const char* const name_;
    const int line_;
};

#define DEFINE_GRID_KERNEL(name) \
    static constexpr ::GridFunctorMeta functor_description = {__FILE__, __LINE__, name};
#else
#define DEFINE_GRID_KERNEL(name) /* nothing */
#endif

#endif //MICROHHC_CUDA_TILING_H
