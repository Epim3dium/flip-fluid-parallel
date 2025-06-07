#include "particle.hpp"
#include "geometry_func.hpp"
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

void accelerate(Particles& particles, vec2f gravity) {
    for(int i = 0; i < max_particle_count; i++) {
        particles.acceleration[i] += gravity;
    }
}
void integrate(Particles& particles, float dt) {
    for(int i = 0; i < max_particle_count; i++) {
        particles.velocity[i] += particles.acceleration[i] * dt;
        particles.position[i] += particles.velocity[i] * dt;
        particles.acceleration[i] = {0, 0};
    }
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
    for(int i = 0; i < max_particle_count; i++) {
        for(int ii = i+1; ii < max_particle_count; ii++) {
        }
    }
}
struct CompactVec {
    uint32_t data[32];
    uint16_t size = 0;
};

__global__ void compareWithNeighbours(Particles& particles, uint32_t* active, uint32_t active_size, int max_segs_cols, const CompactVec* grid) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= active_size)
        return;
    auto row = active[i] / max_segs_cols;
    auto col = active[i] % max_segs_cols;
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
void collide(Particles& particles, AABB sim_area) {
    const uint32_t max_segs_cols = sim_area.size().x / particles.diameter + 1 + 2;
    const uint32_t max_segs_rows = sim_area.size().y / particles.diameter + 1 + 2;
    static int col_grid_size;
    static CompactVec* col_grid = nullptr;

    if(col_grid_size != max_segs_rows * max_segs_cols) {
        if(col_grid_size != 0)
            cudaFree(col_grid);
        CUDA_CALL(cudaMallocManaged(&col_grid, sizeof(CompactVec) * max_segs_rows * max_segs_cols));
    }
    std::unordered_set<uint32_t> active_containers[4];
    int counter = 0;
    for(int i = 0; i < max_particle_count; i++) {
        uint32_t col = (particles.position[i].x - sim_area.min.x) / particles.diameter;
        uint32_t row = (particles.position[i].y - sim_area.min.y) / particles.diameter;
        if(col+1 >= max_segs_cols || row+1 >= max_segs_rows) {
            continue;
        }
        auto& comp_vec = col_grid[(row+1) * max_segs_cols + col+1];
        if(comp_vec.size == 32U) continue;
        comp_vec.data[comp_vec.size++] = i;
        active_containers[(row % 2)*2 + col%2].insert((row+1) * max_segs_cols + col+1);
    }

    CUDA_CALL(cudaMemcpy(particles.gpu_velocity, particles.velocity, sizeof(vec2f) * max_particle_count, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(particles.gpu_position, particles.position, sizeof(vec2f) * max_particle_count, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(particles.gpu_acceleration, particles.acceleration, sizeof(vec2f) * max_particle_count, cudaMemcpyHostToDevice));
    size_t max_size = 0;
    for(const auto& v : active_containers)
        max_size = std::max(v.size(), max_size);

    uint32_t* active_idxs;
    CUDA_CALL(cudaMalloc(&active_idxs, sizeof(uint32_t) * max_size));
    for(int i = 0; i < 4; i++) {
        std::vector<uint32_t> active(active_containers[i].begin(),active_containers[i].end());
        CUDA_CALL(cudaMemcpy(active_idxs, active.data(), sizeof(uint32_t) * active.size(), cudaMemcpyHostToDevice));
        auto N = active.size();
        compareWithNeighbours<<<(N + 255) / 256, 256>>>(particles, active_idxs, N, max_segs_cols, col_grid);
        cudaDeviceSynchronize();
    }
    CUDA_CALL(cudaFree(active_idxs));
    CUDA_CALL(cudaMemcpy(particles.velocity, particles.gpu_velocity, sizeof(vec2f) * max_particle_count, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(particles.position, particles.gpu_position, sizeof(vec2f) * max_particle_count, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(particles.acceleration, particles.gpu_acceleration, sizeof(vec2f) * max_particle_count, cudaMemcpyDeviceToHost));

    for(int i = 0; i < 4; i++) {
        for(auto container : active_containers[i]) {
            auto row = container / max_segs_cols;
            auto col = container % max_segs_cols;
            col_grid[row * max_segs_cols + col].size = 0;
        }
    }

}
void constraint(Particles& particles, AABB area) {
    for(int i = 0; i < max_particle_count; i++) {
        if(!isOverlappingPointAABB(particles.position[i], area)) {
            if(particles.position[i].x > area.max.x || particles.position[i].x < area.min.x)
                particles.velocity[i].x = 0;
            if(particles.position[i].y > area.max.y || particles.position[i].y < area.min.y)
                particles.velocity[i].y = 0;
            particles.position[i].x = std::clamp(particles.position[i].x, area.min.x, area.max.x);
            particles.position[i].y = std::clamp(particles.position[i].y, area.min.y, area.max.y);
        }
    }
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

