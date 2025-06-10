#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/Audio.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <numeric>
#include "fluid.hpp"
#include "particle.hpp"
#include "geometry_func.hpp"
#include "time.hpp"


#include <SFML/Window/Mouse.hpp>
#include <iostream>
#include <stdexcept>
using namespace sf;

int main(int argc, char** argv) {
    float w = 1280;
    float h = 720;
    int numParticleIters = 16;
    int numFluidIters = 32;
    float overrelaxation = 1.9f;
    float cell_size_scale = 3.f;
    float raporting_interval = 0.125f;
    bool pushOut = true;
    Particles particles;
    std::string help_msg = R"""(
Usage: fluid-sim [OPTIONS]
Options:
        -w [width] 
        -h [height] 
        -c [cell scale] 
        -r [particle radius] 
        -n [amount of particles]
        -i [num of iterations for particle solver]
        -f [num of iterations for fluid solver]
        -x [overrelaxation coef]
        -d [1/0 should drift be corrected]
        -t [raporting time in seconds]
)""";
    try {
        for(int i = 1; i < argc; i += 2) {
            std::string flag = argv[i];
            std::string arg = argv[i+1];
            if(flag == "-w") {
                w = stoi(arg);
            }else if(flag == "-h") {
                h = stoi(arg);
            }else if(flag == "-c") {
                cell_size_scale = stof(arg);
            }else if(flag == "-r") {
                particles.radius = stof(arg);
            }else if(flag == "-n") {
                particles.max_particle_count = stoi(arg);
            }else if(flag == "-i") {
                numParticleIters = stoi(arg);
            }else if(flag == "-f") {
                numFluidIters = stoi(arg);
            }else if(flag == "-x") {
                overrelaxation = stof(arg);
            }else if(flag == "-d") {
                pushOut = stoi(arg);
            }else if(flag == "-t") {
                raporting_interval = stof(arg);
            }else if(flag == "--help") {
                printf("%s", help_msg.c_str());
                return 0;
            }else {
                throw std::invalid_argument("unrecognized flag");
            }
        }
    }catch(...) {
        printf("incorrect arguments were given!");
        printf("%s", help_msg.c_str());
        return 0;
    }
    auto fluid_cell_size = particles.diameter * cell_size_scale;
    AABB screen_area = AABB::CreateMinSize({0, 0}, {w, h});
    RenderWindow window(VideoMode(w, h), "demo");
    auto area = screen_area;
    area.setSize(area.size() * 0.75f);
    init(particles, area, 1.5f);

    auto fluid_size = screen_area.size() / fluid_cell_size;
    Fluid fluid(fluid_cell_size, fluid_size.x + 1, fluid_size.y + 1);


    std::cout << "{\n";
    auto dispNameValue = [&](std::string name, auto value, bool isLast = false) {
        std::cout << "\t\"" << name << "\" : \"" << value << "\"";
        if(!isLast)
            std::cout << ",";
        std::cout << "\n";
    };
    dispNameValue("particle radius", particles.radius);
    dispNameValue("particle count", Particles::max_particle_count);
    dispNameValue("num of particle iters", numParticleIters);
    dispNameValue("num of fluid iters", numFluidIters);
    dispNameValue("overrelaxation", overrelaxation);
    std::cout << "\t\"measurements\" : [\n";

    float total_time = 0;
    Clock deltaClock;
    bool shouldReport = true;
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
        auto posi = sf::Mouse().getPosition(window);
        vec2f mouse_pos  = {(float)posi.x, (float)posi.y};
        mouse_pos.y = window.getSize().y - mouse_pos.y;
        vec2f mouse_dir = mouse_pos - last_mouse_pos;
        const float brush_size = 50.f;
        if(qlen(mouse_dir) != 0 && sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)) {
            for(int i = 0; i < Particles::max_particle_count; i++) {
                auto scalar = length(mouse_dir) * 100.f;
                if(length(mouse_pos - particles.position[i]) < brush_size && scalar > 0.f) {
                    particles.velocity[i] = normal(mouse_dir) * std::clamp(scalar, 0.f, 1000.f);
                }
            }
        }

        float deltaTime = deltaClock.restart().asSeconds();
        total_time += deltaTime;

        static std::vector<std::map<std::string, float>> times;
        auto cur_times = fluid.simulate(particles, screen_area, deltaTime, vec2f(0, -1000.f), numFluidIters, numParticleIters, overrelaxation, pushOut);
        times.push_back(cur_times);
        if(report_clock.getElapsedTime() > raporting_interval && shouldReport) {
            static bool displayed = false;
            std::cout << "\t\t";
            if(displayed) {
                std::cout << ",";
            }
            displayed = true;
            std::cout << "{\n";
            report_clock.restart();
            auto displayAvg = [&](std::string name, auto& vec_table) {
                float sum = 0.f;
                for(auto& tab : vec_table) {
                    sum += tab[name];
                }
                auto avg = sum / vec_table.size();
                std::cout << "\"" << name << "\" : \"" << avg << "\"";
            };
            for(auto [name, val] : times.front()) {
                std::cout << "\t\t\t";
                displayAvg(name, times);
                std::cout << ",\n";
            }
            float sum = 0.f;
            for(const auto& v : times) {
                for(auto [name, val] : v)
                    sum += val;
            }
            float avg = sum / times.size();
            std::cout << "\t\t\t\"FPS\" : \"" << 1.f/avg << "\"";
            std::cout << ",\n";

            std::cout << "\t\t}\n";
            times.clear();
        }

        window.clear();
        static bool drawParticles = false;
        static bool drawGrid= true;
        static bool pressed = 0;
        if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::P)) {
            if(!pressed) {
                drawParticles = !drawParticles ;
            }
            pressed = true;
        }else if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::G)) {
            if(!pressed) {
                drawGrid = !drawGrid;
            }
            pressed = true;
        }else {
            pressed = false;
        }

        if(drawGrid)
            fluid.draw(screen_area, window);
        if(drawParticles)
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
    cleanup(particles);
    std::cout << "\t]\n}\n";

    return 0;
}
