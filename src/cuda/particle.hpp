#ifndef EMP_PARTICLE_HPP
#define EMP_PARTICLE_HPP
#include "vec2.hpp"
#include "AABB.hpp"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <SFML/Graphics/RenderTarget.hpp>
#include <cstdint>
static constexpr uint32_t max_particle_count = 5000;
struct Particles {
    static constexpr float radius = 3.0f;
    static constexpr float diameter = radius*2;
    vec2f *position;
    vec2f *velocity;
    vec2f *acceleration;

    vec2f *gpu_position;
    vec2f *gpu_velocity;
    vec2f *gpu_acceleration;
};
struct ParticleSolveBlock {
    Particles& particles;
    ParticleSolveBlock(Particles& p);
    ~ParticleSolveBlock();
};
void derive(Particles& particles, float dt);
void accelerate(Particles& particles, vec2f gravity);
void integrate(Particles& particles, float dt);
void collide(Particles& particles);
void collide(Particles& particles, AABB sim_area);
void constraint(Particles& particles, AABB area);
void draw(Particles& particles, sf::RenderTarget& window, sf::Color color = sf::Color(255, 255, 255, 255));
void init(Particles& particles, AABB screen_area, float spacing = 1.f, int seed = 42);
void cleanup(Particles& particles);
#endif
