#ifndef EMP_PARTICLE_HPP
#define EMP_PARTICLE_HPP
#include "circle.hpp"
#include <SFML/Graphics/RenderTarget.hpp>
#include <cstdint>
#include <map>
#include <string>
static constexpr uint32_t max_particle_count = 5000;
struct Particles {
    static constexpr float radius = 3.f;
    static constexpr float diameter = radius*2;
    vec2f position[max_particle_count];
    vec2f velocity[max_particle_count];
    vec2f acceleration[max_particle_count];
};
struct ParticleSolveBlock {
    Particles& particles;
    ParticleSolveBlock(Particles& p) : particles(p) {}
    ~ParticleSolveBlock() {}
};
void derive(Particles& particles, float dt);
void accelerate(Particles& particles, vec2f gravity);
void integrate(Particles& particles, float dt);
void collide(Particles& particles);
std::map<std::string, float> collide(Particles& particles, AABB sim_area);
void constraint(Particles& particles, AABB area);
void draw(Particles& particles, sf::RenderTarget& window, sf::Color color = sf::Color(255, 255, 255, 255));
void init(Particles& particles, AABB screen_area, float spacing = 1.f, int seed = 42);
void cleanup(Particles& particles);
#endif
