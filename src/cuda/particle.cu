#include "particle.hpp"
#include "geometry_func.hpp"
#include "time.hpp"
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
#include <execution>
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

__global__ void accelerateKernel(Particles& particles, vec2f gravity, int max_count) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= max_count) return;
    particles.gpu_acceleration[i] += gravity;
}

void accelerate(Particles& particles, vec2f gravity) {
    int threadsPerBlock = 1024;
    int blocks = (max_particle_count + threadsPerBlock - 1) / threadsPerBlock;
    accelerateKernel<<<blocks, threadsPerBlock>>>(particles, gravity, max_particle_count);
    cudaDeviceSynchronize();
    CUDA_KERNEL_CHECK();
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
    CUDA_KERNEL_CHECK();
}
__device__ void resolveVelocities(const vec2f* velocity, vec2f* velocity_dt, int p1, int p2, vec2f normal) {
    auto rel_vel = velocity[p1] - velocity[p2];
    auto rel_vel_normal = dot(rel_vel, normal);
    if (rel_vel_normal > 0) return;
    float restitution = 0.1f;
    float impulse = -(1 + restitution) * rel_vel_normal * 0.5;

        atomicAdd(&velocity_dt[p1].x, normal.x * impulse);
        atomicAdd(&velocity_dt[p1].y, normal.y * impulse);
        atomicAdd(&velocity_dt[p2].x, -normal.x * impulse);
        atomicAdd(&velocity_dt[p2].y, -normal.y * impulse);
    // velocity_dt[p1] += impulse * normal;
    // velocity_dt[p2] -= impulse * normal;
}
__device__ bool isColliding(const vec2f& v1, const vec2f& v2, float radius) {
    float2 d = make_float2(
        v1.x - v2.x,
        v1.y - v2.y
    );
    float distSquared = d.x * d.x + d.y * d.y;
    float combinedRadius = radius * 2.f;
    return distSquared < (combinedRadius * combinedRadius);
}
__device__ void processCollision(const vec2f* position, const vec2f* velocity,vec2f* position_dt, vec2f* velocity_dt , int i, int ii, float radius) {
    auto diff = position[ii] - position[i];
    const float min_dist = radius * 2;
    if(isColliding(position[i], position[ii], radius)) {
        auto l = length(diff);
        auto n = diff / l;
        static const float damping = 0.7f;
        auto c = (min_dist - l) * 0.5f * damping;
        atomicAdd(&position_dt[i].x, -n.x * c);
        atomicAdd(&position_dt[i].y, -n.y * c);
        atomicAdd(&position_dt[ii].x, n.x * c);
        atomicAdd(&position_dt[ii].y, n.y * c);
        resolveVelocities(velocity, velocity_dt, i, ii, -n);
    }
}
void collide(Particles& particles) {
    assert(false);
}
#define COMP_SIZE 128
struct CompactVec {
    uint32_t data[COMP_SIZE];
    int size = 0;
};

__global__ void compareWithNeighbours(
    const vec2f* position, const vec2f* velocity,
    vec2f* position_dt, vec2f* velocity_dt,
    float radius, uint32_t* active, uint32_t* active_size, int max_segs_cols, int max_segs_rows, const CompactVec* grid) {
    int max_size = max_segs_rows * max_segs_cols;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= *active_size)
        return;
    auto row = active[i] / max_segs_cols;
    auto col = active[i] % max_segs_cols;
    auto& comp_vec1 = grid[row*max_segs_cols + col];
    if(comp_vec1.size == 0) 
        return;
    #define offset_grid(dirx, diry) &grid[(row+dirx) * max_segs_cols + col+diry]
    const CompactVec* comp_vecs[4] = {offset_grid(1, 0), offset_grid(0, 1), offset_grid(1, 1),offset_grid(1,-1)};

    for(int i = 0; i < comp_vec1.size; i++) {
        auto idx1 = comp_vec1.data[i];
        for(int ii = i + 1; ii < min(comp_vec1.size, COMP_SIZE); ii++) {
            auto idx2 = comp_vec1.data[ii];
            processCollision(position, velocity, position_dt, velocity_dt, idx1, idx2, radius);
        }
        for(int j = 0; j < 4U; j++) {
            auto other = comp_vecs[j];
            for(auto ii = 0; ii < min(other->size, COMP_SIZE); ii++) {
                auto idx2 = (*other).data[ii];
                processCollision(position, velocity, position_dt, velocity_dt, idx1, idx2, radius);
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
    const vec2f* position, float diameter,
    CompactVec* col_grid,
    vec2f sim_area_min,
    int max_segs_cols,
    int max_segs_rows,
    uint32_t* active, // Flattened grid of flags (per cell)
    uint32_t* active_size
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= max_particle_count) return;

    vec2f pos = position[i];
    uint32_t col = (pos.x - sim_area_min.x) / diameter;
    uint32_t row = (pos.y - sim_area_min.y) / diameter;

    if (col + 1 >= max_segs_cols || row + 1 >= max_segs_rows) return;

    int grid_idx = (row + 1) * max_segs_cols + (col + 1);
    CompactVec& comp_vec = col_grid[grid_idx];

    // Atomically get index to write into comp_vec
    uint32_t insert_idx = atomicAdd(&comp_vec.size, 1);
    if (insert_idx <= COMP_SIZE) {
        comp_vec.data[insert_idx] = i;
    }    
    bool old = (insert_idx != 0);  //only the first one insterts
    if(old == 1) return;

    int max_size = max_segs_cols * max_segs_rows;
    auto idx = atomicAdd(active_size, 1);
    active[idx] = grid_idx;
}
__global__ void addDeltas(
    vec2f* v1, vec2f* dt, uint32_t max_size) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= max_size) return;
    v1[i] += dt[i];
}

std::map<std::string, float> collide(Particles& particles, AABB sim_area) {
    std::map<std::string, float> result;
    const uint32_t max_segs_cols = sim_area.size().x / particles.diameter + 1 + 2;
    const uint32_t max_segs_rows = sim_area.size().y / particles.diameter + 1 + 2;

    static int col_grid_size = 0;
    static CompactVec* gpu_col_grid = nullptr;
    static uint32_t* gpu_active_idxs = nullptr;
    static uint32_t* gpu_active_size = nullptr;

    CUDA_CALL(cudaMemset(particles.gpu_position_dt, 0, sizeof(vec2f) * max_particle_count));
    CUDA_CALL(cudaMemset(particles.gpu_velocity_dt, 0, sizeof(vec2f) * max_particle_count));
    Stopwatch stop;
    if(col_grid_size != max_segs_rows * max_segs_cols) {
        if(col_grid_size != 0) {
            cudaFree(gpu_col_grid);
            cudaFree(gpu_active_idxs);
        }else {
            CUDA_CALL(cudaMalloc(&gpu_active_size, sizeof(uint32_t)));
            CUDA_CALL(cudaMemset(gpu_active_size, 0, sizeof(uint32_t)));
        }
        col_grid_size = max_segs_rows * max_segs_cols;
        CUDA_CALL(cudaMalloc(&gpu_active_idxs, sizeof(uint32_t) * col_grid_size));
        CUDA_CALL(cudaMalloc(&gpu_col_grid, sizeof(CompactVec) * col_grid_size));
        CUDA_CALL(cudaMemset(gpu_col_grid, 0, sizeof(CompactVec) * col_grid_size));
    }
    result["particles::collide::allocate"] += stop.restart();


    int threadsPerBlock = 1024;
    int blocks = (max_particle_count + threadsPerBlock - 1) / threadsPerBlock;
    assignParticlesToGrid<<<blocks, threadsPerBlock>>>(
        particles.gpu_position, particles.diameter,
        gpu_col_grid, 
        sim_area.min,
        max_segs_cols, max_segs_rows,
        gpu_active_idxs, gpu_active_size
    );
    cudaDeviceSynchronize();
    result["particles::collide::assign"] += stop.restart();

    compareWithNeighbours<<<(col_grid_size + threadsPerBlock) / threadsPerBlock, threadsPerBlock>>>(
        particles.gpu_position,
        particles.gpu_velocity,
        particles.gpu_position_dt,
        particles.gpu_velocity_dt,
        particles.radius,
        gpu_active_idxs,
        gpu_active_size,
        max_segs_cols,
        max_segs_rows,
        gpu_col_grid);
    cudaDeviceSynchronize();
    result["particles::collide::compare"] += stop.restart();
    blocks = (max_particle_count + threadsPerBlock - 1) / threadsPerBlock;
    addDeltas<<<blocks, threadsPerBlock>>>(particles.gpu_position, particles.gpu_position_dt, max_particle_count);
    addDeltas<<<blocks, threadsPerBlock>>>(particles.gpu_velocity, particles.gpu_velocity_dt, max_particle_count);
    cudaDeviceSynchronize();

    blocks = (col_grid_size + threadsPerBlock - 1) / threadsPerBlock;
    // clearGrid<<<blocks, threadsPerBlock>>>(gpu_col_grid, col_grid_size);
    CUDA_CALL(cudaMemset(gpu_active_size, 0, sizeof(uint32_t)));
    CUDA_CALL(cudaMemset(gpu_col_grid, 0, sizeof(CompactVec) * col_grid_size));
    CUDA_CALL(cudaMemset(gpu_active_idxs, 0, sizeof(uint32_t) * col_grid_size));
    result["particles::collide::reset"] += stop.restart();
    return result;
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
    CUDA_KERNEL_CHECK();
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

    CUDA_CALL(cudaMalloc(&particles.gpu_velocity_dt, sizeof(vec2f) * max_particle_count));
    CUDA_CALL(cudaMalloc(&particles.gpu_position_dt, sizeof(vec2f) * max_particle_count));

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

