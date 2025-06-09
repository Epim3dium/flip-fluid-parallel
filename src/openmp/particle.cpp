#include "particle.hpp"
#include "geometry_func.hpp"
#include "time.hpp"
#include <SFML/Graphics/CircleShape.hpp>
#include <array>
#include <cmath>
#include <cstdlib>
#include <execution>
#include <iostream>
#include <stdexcept>
#include <unordered_set>
#include <omp.h>

void cleanup(Particles& particles) {
}
void accelerate(Particles& particles, vec2f gravity) {
    #pragma omp parallel for
    for(int i = 0; i < max_particle_count; i++) {
        particles.acceleration[i] += gravity;
    }
}
void integrate(Particles& particles, float dt) {
    #pragma omp parallel for
    for(int i = 0; i < max_particle_count; i++) {
        particles.velocity[i] += particles.acceleration[i] * dt;
        particles.position[i] += particles.velocity[i] * dt;
        particles.acceleration[i] = {0, 0};
    }
}
void resolveVelocities(Particles& particles, int p1, int p2, vec2f normal) {
    auto rel_vel = particles.velocity[p1] - particles.velocity[p2];
    auto rel_vel_normal = dot(rel_vel, normal);
    if (rel_vel_normal > 0) return;
    float restitution = 0.1f;
    float impulse = -(1 + restitution) * rel_vel_normal * 0.5;

    particles.velocity[p1] += impulse * normal;
    particles.velocity[p2] -= impulse * normal;

}
void processCollision(Particles& particles, int i, int ii) {
    auto diff = particles.position[ii] - particles.position[i];
    const float min_dist = particles.radius * 2;
    auto ql = qlen(diff);
    if(ql < min_dist * min_dist && ql > 1e-10f) {
        auto l = length(diff);
        auto n = diff / l;
        static const float damping = 0.7f;
        auto c = (min_dist - l) * 0.5f * damping;
        particles.position[i] -= n * c;
        particles.position[ii] += n * c;
        // resolveVelocities(particles, i, ii, -n);
    }
}
void collide(Particles& particles) {
    for(int i = 0; i < max_particle_count; i++) {
        for(int ii = i+1; ii < max_particle_count; ii++) {
            processCollision(particles, i, ii);
        }
    }
}
typedef std::vector<uint32_t> CompactVec;

void compareWithNeighbours(Particles& particles, int col, int row, int max_segs_rows, const std::vector<CompactVec>& grid) {
    auto& comp_vec1 = grid[row*max_segs_rows + col];
    if(comp_vec1.size() == 0) return;
    #define offset_grid(dirx, diry) &grid[(row+dirx) * max_segs_rows + col+diry]
    std::array<const CompactVec*, 4U> comp_vecs = {offset_grid(1, 0), offset_grid(0, 1), offset_grid(1, 1),offset_grid(1,-1)};

    for(int i = 0; i < comp_vec1.size(); i++) {
        auto idx1 = comp_vec1[i];
        for(int ii = i + 1; ii < comp_vec1.size(); ii++) {
            auto idx2 = comp_vec1[ii];
            processCollision(particles, idx1, idx2);
        }
        for(auto other : comp_vecs) {
            for(auto ii = 0; ii < other->size(); ii++) {
                auto idx2 = (*other)[ii];
                processCollision(particles, idx1, idx2);
            }
        }
    }
}
std::map<std::string, float> collide(Particles& particles, AABB sim_area) {
    std::map<std::string, float> result;
    const uint32_t max_segs_cols = sim_area.size().x / particles.diameter + 1 + 2;
    const uint32_t max_segs_rows = sim_area.size().y / particles.diameter + 1 + 2;
    static std::vector<CompactVec> col_grid;

    Stopwatch stop;
    if(col_grid.size() != max_segs_rows * max_segs_cols) {
        col_grid = std::vector<CompactVec>(max_segs_rows*max_segs_cols);
    }
    auto max_dim = std::max(max_segs_cols, max_segs_rows);
    std::unordered_set<uint32_t> active_containers[4];
    int counter = 0;
    for(int i = 0; i < max_particle_count; i++) {
        uint32_t col = (particles.position[i].x - sim_area.min.x) / particles.diameter;
        uint32_t row = (particles.position[i].y - sim_area.min.y) / particles.diameter;
        if(col+1 >= max_segs_cols || row+1 >= max_segs_rows) {
            std::cerr << "out of range particle\n";
            continue;
        }
        auto& comp_vec = col_grid[(row+1) * max_segs_rows + col+1];
        comp_vec.push_back(i);
        active_containers[(row % 2)*2 + col%2].insert((row+1) * max_dim + col+1);
    }
    result["particles::collide::assign"] += stop.restart();

    for(int i = 0; i < 4; i++) {
        std::vector<uint32_t> active(active_containers[i].begin(),active_containers[i].end());
        #pragma omp parallel for
        for(int i = 0; i < active.size(); i++) {
            auto row = active[i] / max_dim;
            auto col = active[i] % max_dim;
            compareWithNeighbours(particles, col, row, max_segs_rows, col_grid);
        }
    }
    result["particles::collide::compare"] += stop.restart();

    for(int i = 0; i < 4; i++) {
        for(auto container : active_containers[i]) {
            auto row = container / max_dim;
            auto col = container % max_dim;
            col_grid[row * max_segs_rows + col].clear();
        }
    }
    result["particles::collide::cleanup"] += stop.restart();
    return result;
}
void constraint(Particles& particles, AABB area) {
    #pragma omp parallel for
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
void draw(Particles& particles, sf::RenderTarget& window, sf::Color color) {
    sf::VertexArray quads(sf::Quads, 4 * max_particle_count);

    for(int i = 0; i < max_particle_count; i += 1) { 
        auto pos = particles.position[i];
        pos.y = window.getSize().y - pos.y;
        // define the position of the triangle's points
        quads[i * 4 + 0].position = sf::Vector2f(pos.x, pos.y) + sf::Vector2f(0.f, particles.radius);
        quads[i * 4 + 1].position = sf::Vector2f(pos.x, pos.y) + sf::Vector2f(particles.radius, 0);
        quads[i * 4 + 2].position = sf::Vector2f(pos.x, pos.y) + sf::Vector2f(0.f, -particles.radius);
        quads[i * 4 + 3].position = sf::Vector2f(pos.x, pos.y) + sf::Vector2f(-particles.radius, 0);

        quads[i * 4 + 0].color = color;
        quads[i * 4 + 1].color = color;
        quads[i * 4 + 2].color = color;
        quads[i * 4 + 3].color = color;
    }
    window.draw(quads);
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

