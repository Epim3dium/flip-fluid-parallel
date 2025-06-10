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
void draw(Particles& particles, sf::RenderTarget& window, sf::Color color) {
    sf::VertexArray quads(sf::Quads, 4 * Particles::max_particle_count);

    for(int i = 0; i < Particles::max_particle_count; i += 1) { 
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

