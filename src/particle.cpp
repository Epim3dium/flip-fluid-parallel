#include "particle.hpp"
#include "geometry_func.hpp"
#include <SFML/Graphics/CircleShape.hpp>
#include <cstdlib>
#include <execution>
void derive(Particles& particles, float dt) {
    for(int i = 0; i < max_particle_count; i++) {
        particles.velocity[i] = (particles.position[i] - particles.prev_position[i]);
        particles.prev_position[i] = particles.position[i];
    }
}

void accelerate(Particles& particles, vec2f gravity) {
    for(int i = 0; i < max_particle_count; i++) {
        particles.acceleration[i] += gravity;
    }
}
void integrate(Particles& particles, float dt) {
    for(int i = 0; i < max_particle_count; i++) {
        particles.position[i] += particles.velocity[i] + particles.acceleration[i] * dt * dt;;
        particles.acceleration[i] = {0, 0};
    }
}
void collide(Particles& particles) {
    for(int i = 0; i < max_particle_count; i++) {
        for(int ii = i+1; ii < max_particle_count; ii++) {
            auto diff = particles.position[ii] - particles.position[i];
            auto min_dist = particles.radius[i] + particles.radius[ii];
            if(qlen(diff) < min_dist * min_dist) {
                auto l = length(diff);
                auto c = min_dist - l;
                auto wi = particles.radius[ii] / min_dist;
                auto wii = particles.radius[i] / min_dist;
                particles.position[i] -= normal(diff) * c * wi;
                particles.position[ii] += normal(diff) * c * wii;
            }
        }
    }
}
void constraint(Particles& particles, AABB area) {
    for(int i = 0; i < max_particle_count; i++) {
        particles.position[i].x = std::clamp(particles.position[i].x, area.min.x, area.max.x);
        particles.position[i].y = std::clamp(particles.position[i].y, area.min.y, area.max.y);
    }
}
void draw(Particles& particles, sf::RenderTarget& window) {
    sf::CircleShape circ;
    circ.setPointCount(16U);
    for(int i = 0; i < max_particle_count; i++) { 
        circ.setPosition(particles.position[i].x, window.getSize().y - particles.position[i].y);
        auto rad = particles.radius[i];
        circ.setRadius(rad);
        circ.setOrigin(rad, rad);
        circ.setFillColor(particles.color[i]);
        window.draw(circ);
    }
}
void init_random(Particles& particles, float screen_width, std::pair<float, float> velocity_range, std::pair<float, float> radius_range, int seed) {
    srand(seed);
    screen_width -= radius_range.second;
    int width = screen_width / (radius_range.second * 2);
    for(int i = 0; i < max_particle_count; i++) { 
        char r = ((float)rand() / RAND_MAX) * 255;
        particles.color[i].r = 255;
        particles.color[i].g = 255;
        particles.color[i].b = 255;
        float rad = ((float)rand() / RAND_MAX);
        rad = (rad * (radius_range.second - radius_range.first) + radius_range.first);
        particles.radius[i] = rad;

        float offsetx = ((i/width) % 2) * radius_range.second;
        particles.position[i].x = (i % width) * radius_range.second * 2 + offsetx;
        particles.position[i].y = (i / width) * radius_range.second * 2;
        particles.prev_position[i].x = particles.position[i].x;
        particles.prev_position[i].y = particles.position[i].y;

        float vx= (((float)rand() / RAND_MAX) - 0.5) * 2.f;
        float vy= (((float)rand() / RAND_MAX) - 0.5) * 2.f;
        vec2f n = normal({vx, vy});
        float t = (float)rand() / RAND_MAX;
        particles.velocity[i].x = n.x * ( t * (velocity_range.second - velocity_range.first) + velocity_range.first);
        particles.velocity[i].y = n.y * ( t * (velocity_range.second - velocity_range.first) + velocity_range.first);
        particles.acceleration[i].x = 0;
        particles.acceleration[i].y = 0;
    }
}

