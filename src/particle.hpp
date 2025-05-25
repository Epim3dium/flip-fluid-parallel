#ifndef EMP_PARTICLE_HPP
#define EMP_PARTICLE_HPP
#include "circle.hpp"
#include <SFML/Graphics/RenderTarget.hpp>
#include <cstdint>
static constexpr uint32_t max_particle_count = 1101U;
struct Particles {
    vec2f prev_position[max_particle_count];
    vec2f position[max_particle_count];
    vec2f velocity[max_particle_count];
    vec2f acceleration[max_particle_count];
    float radius[max_particle_count];
    Color color[max_particle_count];
};
void derive(Particles& particles, float dt);
void accelerate(Particles& particles, vec2f gravity);
void integrate(Particles& particles, float dt);
void collide(Particles& particles);
void constraint(Particles& particles, AABB area);
void draw(Particles& particles, sf::RenderTarget& window);
void init_random(Particles& particles, float width, std::pair<float, float> velocity_range, std::pair<float, float> radius_range, int seed);
#endif
