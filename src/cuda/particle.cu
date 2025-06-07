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
    auto rel_vel = particles.velocity[p1] - particles.velocity[p2];
    auto rel_vel_normal = dot(rel_vel, normal);
    if (rel_vel_normal > 0) return;
    float restitution = 0.1f;
    float impulse = -(1 + restitution) * rel_vel_normal * 0.5;

    particles.velocity[p1] += impulse * normal;
    particles.velocity[p2] -= impulse * normal;

}
__device__ void processCollision(Particles& particles, int i, int ii) {
    auto diff = particles.position[ii] - particles.position[i];
    const float min_dist = particles.radius * 2;
    if(qlen(diff) < min_dist * min_dist && qlen(diff) > 1e-10f) {
        auto l = length(diff);
        auto n = normal(diff);
        static const float damping = 0.7f;
        auto c = (min_dist - l) * 0.5f * damping;
        particles.position[i] -= n * c;
        particles.position[ii] += n * c;
        resolveVelocities(particles, i, ii, -n);
    }
}
__device__ bool isColliding(const Particles& particles, int idx1, int idx2) {
    float2 d = make_float2(
        particles.position[idx1].x - particles.position[idx2].x,
        particles.position[idx1].y - particles.position[idx2].y
    );
    float distSquared = d.x * d.x + d.y * d.y;
    float combinedRadius = particles.radius * 2.f;
    return distSquared < (combinedRadius * combinedRadius);
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

__global__ void compareWithNeighbours(Particles& particles, int col, int row, int max_segs_rows, const CompactVec* grid) {
    auto& comp_vec1 = grid[row*max_segs_rows + col];
    if(comp_vec1.size == 0) return;
    #define offset_grid(dirx, diry) &grid[(row+dirx) * max_segs_rows + col+diry]
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
    static CompactVec* col_grid;
    static size_t col_grid_size = 0;

    if(col_grid_size != max_segs_rows * max_segs_cols) {
        if(col_grid_size != 0)
            cudaFree(col_grid);
        CUDA_CALL(cudaMallocManaged(&col_grid, sizeof(CompactVec) * max_segs_rows * max_segs_cols));
    }
    auto max_dim = std::max(max_segs_cols, max_segs_rows);
    std::unordered_set<uint32_t> active_containers;
    for(int i = 0; i < max_particle_count; i++) {
        uint32_t col = (particles.position[i].x - sim_area.min.x) / particles.diameter;
        uint32_t row = (particles.position[i].y - sim_area.min.y) / particles.diameter;
        if(col+1 >= max_segs_cols || row+1 >= max_segs_rows) {
            continue;
        }
        auto& comp_vec = col_grid[(row+1) * max_segs_rows + col+1];
        if(comp_vec.size == 32U) continue;
        comp_vec.data[comp_vec.size++] = i;
        active_containers.insert((row+1) * max_dim + col+1);
    }

    for(auto i : active_containers) {
        auto row = i / max_dim;
        auto col = i % max_dim;
        compareWithNeighbours<<<1, 1>>>(particles, col, row, max_segs_rows, col_grid);
    }

    for(auto container : active_containers) {
        auto row = container / max_dim;
        auto col = container % max_dim;
        col_grid[row * max_segs_rows + col].size = 0;
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
    srand(seed);
    auto w = screen_area.size().x - particles.radius * 2.f;
    int width = w / (particles.radius * 2 * spacing);
    for(int i = 0; i < max_particle_count; i++) { 
        particles.position[i].x = (i % width) * particles.radius * 2.f * spacing + screen_area.min.x;
        particles.position[i].y = (i / width) * particles.radius * 2.f * spacing + screen_area.min.y;

        particles.velocity[i] =  {0, 0};
        particles.acceleration[i] = {0, 0};
    }
}

