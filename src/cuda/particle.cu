#include "particle.hpp"
#include "geometry_func.hpp"
#include <cassert>
#include <device_atomic_functions.h>
#include "vec2.hpp"
#include "AABB.hpp"
#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/RenderStates.hpp>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cuda_device_runtime_api.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>
#include <unordered_set>
#define CUDA_KERNEL_CHECK()                                                     \
do {                                                                            \
    cudaError_t err = cudaGetLastError();                                       \
    if (err != cudaSuccess) {                                                   \
        fprintf(stderr, "Kernel launch error at %s:%d: %s\n",                   \
                __FILE__, __LINE__, cudaGetErrorString(err));                   \
        exit(EXIT_FAILURE);                                                     \
    }                                                                           \
} while (0)
#define CUDA_CALL(call)                                                         \
do {                                                                            \
    cudaError_t err = call;                                                     \
    if (err != cudaSuccess) {                                                   \
        fprintf(stderr, "CUDA error in %s (%s:%d): %s\n",                       \
                #call, __FILE__, __LINE__, cudaGetErrorString(err));            \
        exit(EXIT_FAILURE);                                                     \
    }                                                                           \
} while (0)

__global__ void accelerateKernel(Particles& particles, vec2f gravity) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= max_particle_count) return;
    particles.gpu_acceleration[i] += gravity;
}

void accelerate(Particles& particles, vec2f gravity) {
    int threadsPerBlock = 1024;
    int blocks = (max_particle_count + threadsPerBlock - 1) / threadsPerBlock;
    accelerateKernel<<<blocks, threadsPerBlock>>>(particles, gravity);
    cudaDeviceSynchronize();
}
__global__ void integrateKernel(Particles& particles, float dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= max_particle_count) return;
    particles.gpu_velocity[i] += particles.gpu_acceleration[i] * dt;
    particles.gpu_position[i] += particles.gpu_velocity[i] * dt;
    particles.gpu_acceleration[i] = {0, 0};
}
void integrate(Particles& particles, float dt) {
    int threadsPerBlock = 1024;
    int blocks = (max_particle_count + threadsPerBlock - 1) / threadsPerBlock;
    integrateKernel<<<blocks, threadsPerBlock>>>(particles, dt);
    cudaDeviceSynchronize();
}
__device__ void resolveVelocities(Particles& particles, int p1, int p2, vec2f normal) {
    auto rel_vel = particles.gpu_velocity[p1] - particles.gpu_velocity[p2];
    auto rel_vel_normal = dot(rel_vel, normal);
    if (rel_vel_normal > 0) return;
    float restitution = 0.1f;
    float impulse = -(1 + restitution) * rel_vel_normal * 0.5;

    particles.gpu_velocity[p1] += impulse * normal;
    particles.gpu_velocity[p2] -= impulse * normal;

}
__device__ bool isColliding(const Particles& particles, int idx1, int idx2) {
    float2 d = make_float2(
        particles.gpu_position[idx1].x - particles.gpu_position[idx2].x,
        particles.gpu_position[idx1].y - particles.gpu_position[idx2].y
    );
    float distSquared = d.x * d.x + d.y * d.y;
    float combinedRadius = particles.radius * 2.f;
    return distSquared < (combinedRadius * combinedRadius);
}
__device__ void processCollision(Particles& particles, int i, int ii) {
    auto diff = particles.gpu_position[ii] - particles.gpu_position[i];
    const float min_dist = particles.radius * 2;
    if(isColliding(particles, i, ii)) {
        auto l = length(diff);
        auto n = normal(diff);
        static const float damping = 0.7f;
        auto c = (min_dist - l) * 0.5f * damping;
        particles.gpu_position[i] -= n * c;
        particles.gpu_position[ii] += n * c;
        resolveVelocities(particles, i, ii, -n);
    }
}
//launch only one block
__global__ void detectCollisionKernel(Particles& particles, bool* collisions) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= max_particle_count) return;

    for (int j = 0; j < max_particle_count; ++j) {
        if (i != j && isColliding(particles, i, j)) {
            collisions[i * max_particle_count + j] = true; // Store pairwise collision result
        }
    }
}
void collide(Particles& particles) {
    assert(false);
}
struct CompactVec {
    uint32_t data[32];
    uint32_t size = 0;
};

__global__ void compareWithNeighbours(Particles& particles, uint32_t* active, uint32_t active_size, uint32_t checkerboard, int max_segs_cols, int max_segs_rows, const CompactVec* grid) {
    int max_size = max_segs_rows * max_segs_cols;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= active_size)
        return;
    auto row = active[i + max_size * checkerboard] / max_segs_cols;
    auto col = active[i + max_size * checkerboard] % max_segs_cols;
    auto& comp_vec1 = grid[row*max_segs_cols + col];
    if(comp_vec1.size == 0) return;
    #define offset_grid(dirx, diry) &grid[(row+dirx) * max_segs_cols + col+diry]
    const CompactVec* comp_vecs[4] = {offset_grid(1, 0), offset_grid(0, 1), offset_grid(1, 1),offset_grid(1,-1)};

    for(int i = 0; i < comp_vec1.size; i++) {
        auto idx1 = comp_vec1.data[i];
        for(int ii = i + 1; ii < comp_vec1.size; ii++) {
            auto idx2 = comp_vec1.data[ii];
            processCollision(particles, idx1, idx2);
        }
        for(auto other : comp_vecs) {
            for(auto ii = 0; ii < other->size; ii++) {
                auto idx2 = (*other).data[ii];
                processCollision(particles, idx1, idx2);
            }
        }
    }
}
__global__ void clearGrid(
    CompactVec* col_grid,
    int col_grid_size) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= col_grid_size) return;
    col_grid[i].size = 0;
}
__global__ void assignParticlesToGrid(
    Particles& particles,
    CompactVec* col_grid,
    int* grid_flags,
    vec2f sim_area_min,
    int max_segs_cols,
    int max_segs_rows,
    uint32_t* active_flags, // Flattened grid of flags (per cell)
    uint32_t active_flags_sizes[4]
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= max_particle_count) return;

    vec2f pos = particles.gpu_position[i];
    uint32_t col = (pos.x - sim_area_min.x) / particles.diameter;
    uint32_t row = (pos.y - sim_area_min.y) / particles.diameter;

    if (col + 1 >= max_segs_cols || row + 1 >= max_segs_rows) return;

    int grid_idx = (row + 1) * max_segs_cols + (col + 1);
    CompactVec& comp_vec = col_grid[grid_idx];

    // Atomically get index to write into comp_vec
    uint32_t insert_idx = atomicAdd(&comp_vec.size, 1);
    if (insert_idx < 32) {
        comp_vec.data[insert_idx] = i;
    }
    bool old = atomicExch(&grid_flags[grid_idx], 1);  // atomically set to 1
    if(old == 1) return;

    int checkerboard = (row % 2)*2 + col%2;
    int max_size = max_segs_cols * max_segs_rows;
    auto idx = atomicAdd(&active_flags_sizes[checkerboard], 1);
    active_flags[idx+checkerboard*max_size] = grid_idx;
}
void collide(Particles& particles, AABB sim_area) {
    const uint32_t max_segs_cols = sim_area.size().x / particles.diameter + 1 + 2;
    const uint32_t max_segs_rows = sim_area.size().y / particles.diameter + 1 + 2;

    static int col_grid_size = 0;
    static CompactVec* gpu_col_grid = nullptr;
    static int* gpu_grid_flags = nullptr;
    static uint32_t* gpu_active_idxs = nullptr;

    if(col_grid_size != max_segs_rows * max_segs_cols) {
        if(col_grid_size != 0) {
            cudaFree(gpu_col_grid);
            cudaFree(gpu_active_idxs);
            cudaFree(gpu_grid_flags);
        }
        col_grid_size = max_segs_rows * max_segs_cols;
        CUDA_CALL(cudaMalloc(&gpu_active_idxs, sizeof(uint32_t) * col_grid_size * 4U));
        CUDA_CALL(cudaMalloc(&gpu_col_grid, sizeof(CompactVec) * col_grid_size));
        CUDA_CALL(cudaMemset(gpu_col_grid, 0, sizeof(CompactVec) * col_grid_size));
        CUDA_CALL(cudaMalloc(&gpu_grid_flags, sizeof(int) * col_grid_size));
        CUDA_CALL(cudaMemset(gpu_grid_flags, 0, sizeof(int) * col_grid_size));
    }

    uint32_t active_flags_sizes[4] = {0, 0, 0, 0};
    int threadsPerBlock = 1024;
    int blocks = (max_particle_count + threadsPerBlock - 1) / threadsPerBlock;
    assignParticlesToGrid<<<blocks, threadsPerBlock>>>(
        particles, gpu_col_grid, gpu_grid_flags,
        sim_area.min,
        max_segs_cols, max_segs_rows,
        gpu_active_idxs, active_flags_sizes
    );
    cudaDeviceSynchronize();

    for(int i = 0; i < 4; i++) {
        auto N = active_flags_sizes[i];
        compareWithNeighbours<<<(N + 1024) / 1024, 1024>>>(particles, gpu_active_idxs, N, i, max_segs_cols, max_segs_rows, gpu_col_grid);
        cudaDeviceSynchronize();
    }

    blocks = (col_grid_size + threadsPerBlock - 1) / threadsPerBlock;
    clearGrid<<<blocks, threadsPerBlock>>>(gpu_col_grid, col_grid_size);
    CUDA_CALL(cudaMemset(gpu_grid_flags, 0, sizeof(int) * col_grid_size));
}
__global__ void constraintKernel(Particles& particles, vec2f area_min, vec2f area_max) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= max_particle_count) return;
    if(particles.gpu_position[i].x > area_max.x || particles.gpu_position[i].x < area_min.x)
        particles.gpu_velocity[i].x = 0;
    if(particles.gpu_position[i].y > area_max.y || particles.gpu_position[i].y < area_min.y)
        particles.gpu_velocity[i].y = 0;
    particles.gpu_position[i].x = fmaxf(particles.gpu_position[i].x, area_min.x);
    particles.gpu_position[i].y = fmaxf(particles.gpu_position[i].y, area_min.y);
    particles.gpu_position[i].x = fminf(particles.gpu_position[i].x, area_max.x);
    particles.gpu_position[i].y = fminf(particles.gpu_position[i].y, area_max.y);
}
void constraint(Particles& particles, AABB area) {
    int threadsPerBlock = 1024;
    int blocks = (max_particle_count + threadsPerBlock - 1) / threadsPerBlock;
    constraintKernel<<<blocks, threadsPerBlock>>>(particles, area.min, area.max);
    cudaDeviceSynchronize();
}
ParticleSolveBlock::ParticleSolveBlock(Particles& p) : particles(p) {
    CUDA_CALL(cudaMemcpy(particles.gpu_velocity, particles.velocity, sizeof(vec2f) * max_particle_count, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(particles.gpu_position, particles.position, sizeof(vec2f) * max_particle_count, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(particles.gpu_acceleration, particles.acceleration, sizeof(vec2f) * max_particle_count, cudaMemcpyHostToDevice));
}
ParticleSolveBlock::~ParticleSolveBlock() {
    CUDA_CALL(cudaMemcpy(particles.velocity, particles.gpu_velocity, sizeof(vec2f) * max_particle_count, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(particles.position, particles.gpu_position, sizeof(vec2f) * max_particle_count, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(particles.acceleration, particles.gpu_acceleration, sizeof(vec2f) * max_particle_count, cudaMemcpyDeviceToHost));
}
void init(Particles& particles, AABB screen_area, float spacing, int seed) {
    auto w = screen_area.size().x - particles.radius * 2.f;
    int width = w / (particles.radius * 2 * spacing);
    CUDA_CALL(cudaMallocHost(&particles.position, sizeof(vec2f) * max_particle_count));
    CUDA_CALL(cudaMalloc(&particles.gpu_position, sizeof(vec2f) * max_particle_count));
    CUDA_CALL(cudaMallocHost(&particles.acceleration, sizeof(vec2f) * max_particle_count));
    CUDA_CALL(cudaMalloc(&particles.gpu_acceleration, sizeof(vec2f) * max_particle_count));
    CUDA_CALL(cudaMallocHost(&particles.velocity, sizeof(vec2f) * max_particle_count));
    CUDA_CALL(cudaMalloc(&particles.gpu_velocity, sizeof(vec2f) * max_particle_count));

    for(int i = 0; i < max_particle_count; i++) { 
        particles.position[i].x = (i % width) * particles.radius * 2.f * spacing + screen_area.min.x;
        particles.position[i].y = (i / width) * particles.radius * 2.f * spacing + screen_area.min.y;

        particles.velocity[i] =  {0, 0};
        particles.acceleration[i] = {0, 0};
    }
}
void cleanup(Particles& particles) {
    CUDA_CALL(cudaFreeHost(particles.position));
    CUDA_CALL(cudaFree(particles.gpu_position));
    CUDA_CALL(cudaFreeHost(particles.acceleration));
    CUDA_CALL(cudaFree(particles.gpu_acceleration));
    CUDA_CALL(cudaFreeHost(particles.velocity));
    CUDA_CALL(cudaFree(particles.gpu_velocity));
}

