#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/Audio.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include "particle.hpp"
#include "time.hpp"


#include <SFML/Window/Mouse.hpp>
#include <iostream>
using namespace sf;

int main() {
    float w = 1280;
    float h = 720;
    AABB screen_area = AABB::CreateMinSize({0, 0}, {w, h});
    RenderWindow window(VideoMode(w, h), "demo");
    window.setFramerateLimit(60);
    Particles particles;
    init_random(particles, w, {0, 0}, {7, 10}, 42);

    float total_time = 0;
    Clock deltaClock;
    Stopwatch report_clock;
    report_clock.restart();
    while (window.isOpen()) {
        Event event;
        while (window.pollEvent(event)) {

            if (event.type == Event::Closed) {
                window.close();
            }
        }
        static vec2f last_mouse_pos;
        vec2f mouse_pos  = (vec2f)sf::Mouse().getPosition(window);
        mouse_pos.y = window.getSize().y - mouse_pos.y;
        vec2f mouse_dir = mouse_pos - last_mouse_pos;
        if(qlen(mouse_dir) != 0) {
            for(int i = 0; i < max_particle_count; i++) {
                if(length(mouse_pos - particles.position[i]) < 10.f) {
                    particles.position[i] += normal(mouse_dir) * qlen(mouse_dir) * 0.01f;
                }
            }
        }

        float deltaTime = deltaClock.restart().asSeconds();
        total_time += deltaTime;

        constexpr int steps = 8;
        float dt = deltaTime / (float)steps;
        accelerate(particles,  vec2f(0, -6000.f));
        float solver_time = 0.f;
        float collision_time = 0.f;
        for(int i = 0; i < steps; i++) {
            {
                Stopwatch counter;
                derive(particles, dt);
                integrate(particles, dt);
                solver_time += counter.getElapsedTime();
            }
            {
                Stopwatch counter;
                collide(particles);
                constraint(particles, screen_area);
                collision_time += counter.getElapsedTime();
            }
        }
        if(report_clock.getElapsedTime() > 1.f) {
            report_clock.restart();
            std::cout << "Raport at: " << total_time << "\n";
            std::cout << "\tsolve time: " << solver_time / (solver_time + collision_time) * 100 << "%\n";
            std::cout << "\tsolve time real: " << solver_time * 1000 << "ms\n";
            std::cout << "\tcollision time: " << collision_time / (solver_time + collision_time) * 100 << "%\n";
            std::cout << "\tcollision time real: " << collision_time * 1000 << "ms\n";
        }

        window.clear();
        draw(particles, window);
        window.display();
        last_mouse_pos = mouse_pos;
    }

    return 0;
}
