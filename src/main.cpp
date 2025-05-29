#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/Audio.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <numeric>
#include <omp.h>
#include "fluid.hpp"
#include "particle.hpp"
#include "geometry_func.hpp"
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
    auto fluid_cell_size = particles.diameter * 2.f;
    auto fluid_size = sf::Vector2<int>(screen_area.size() / fluid_cell_size);
    Fluid fluid(fluid_cell_size, fluid_size.x + 1, fluid_size.y + 1);

    int numParticleIters = 16;
    int numFluidIters = 16;
    float overrelaxation = 1.9f;
    auto area = screen_area;
    area.setSize(area.size() * 0.9f);
    init_random(particles, area, 1.1f);

    std::cout << "{\n";
    auto dispNameValue = [&](std::string name, auto value, bool isLast = false) {
        std::cout << "\t\"" << name << "\" : \"" << value << "\"";
        if(!isLast)
            std::cout << ",";
        std::cout << "\n";
    };
    dispNameValue("particle radius", particles.radius);
    dispNameValue("particle count", max_particle_count);
    dispNameValue("num of particle iters", numParticleIters);
    dispNameValue("num of fluid iters", numFluidIters);
    dispNameValue("overrelaxation", overrelaxation, true);
    std::cout << "\t\"measurements\" : [\n";

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
        const float brush_size = 50.f;
        if(qlen(mouse_dir) != 0 && sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)) {
            for(int i = 0; i < max_particle_count; i++) {
                auto scalar = length(mouse_dir) * 100.f;
                if(length(mouse_pos - particles.position[i]) < brush_size && scalar > 0.f) {
                    particles.velocity[i] = normal(mouse_dir) * std::clamp(scalar, 0.f, 1000.f);
                }
            }
        }

        float deltaTime = deltaClock.restart().asSeconds();
        total_time += deltaTime;

        static std::vector<float> fpss;
        static std::vector<float> col_times;
        static std::vector<float> fluid_times;
        auto [col, fl] = fluid.simulate(particles, screen_area, deltaTime, vec2f(0, -1000.f), numFluidIters, numParticleIters, overrelaxation, true);
        col_times.push_back(col);
        fluid_times.push_back(fl);
        fpss.push_back(1.f / deltaTime);
        if(report_clock.getElapsedTime() > 1.f) {
            static bool displayed = false;
            std::cout << "\t\t";
            if(displayed) {
                std::cout << ",";
            }
            displayed = true;
            std::cout << "{\n";
            report_clock.restart();
            auto displayAvg = [&](std::string name, auto& array, bool isLast = false) {
                auto avg = std::reduce(array.begin(), array.end()) / (float)array.size();
                std::cout << "\t\t\t\"" << name << "\" : \"" << avg << "\"";
                if(!isLast)
                    std::cout << ",";
                std::cout << "\n";
            };
            displayAvg("FPS", fpss);
            displayAvg("particle time", col_times);
            displayAvg("fluid time", fluid_times, true);
            std::cout << "\t\t}\n";
            fpss.clear();
            col_times.clear();
            fluid_times.clear();
        }

        window.clear();
        fluid.draw(screen_area, window);
        draw(particles, window, Color(70, 70, 250));
        sf::CircleShape cs(brush_size);
        cs.setOrigin({brush_size, brush_size});
        cs.setPosition(mouse_pos.x, screen_area.size().y - mouse_pos.y);
        cs.setFillColor(Color(0, 0, 0, 0));
        cs.setOutlineColor(Color(255, 255, 255));
        cs.setOutlineThickness(2.f);
        window.draw(cs);
        window.display();
        last_mouse_pos = mouse_pos;
    }
    std::cout << "]\n}\n";

    return 0;
}
