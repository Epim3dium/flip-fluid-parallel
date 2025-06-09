#ifndef EMP_FLUID_HPP
#define EMP_FLUID_HPP
#include "SFML/Graphics/Image.hpp"
#include "SFML/Graphics/Sprite.hpp"
#include "SFML/Graphics/Texture.hpp"
#include "SFML/Graphics/RenderTarget.hpp"
#include "particle.hpp"
#include <algorithm>
#include <iostream>
#include <map>
#include <unordered_map>
enum class eCellTypes {
    Fluid,
    Air,
    Solid
};
class Fluid {
    enum class eFields {
        U,
        V
    };
    int m_width;
    int m_height;
    int m_num_cells;
    float m_cell_size;
    std::vector<Color> cell_color;
    std::vector<eCellTypes> cell_type;
    std::vector<float> solid;
    std::vector<float> pressure;
    std::vector<vec2f> prev_velocities;
    std::vector<vec2f> velocities;
    std::vector<vec2f> velocities_diff;
    std::vector<float> particle_density;
    float particleRestDensity = 0;
public:
    float density = 1;
    float flipRatio = 0.9f;

    void updateParticleDensity(Particles& particles);
    void transferVelocities(bool toGrid, float flipRatio, Particles& particles);
    void solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift = true);

    std::map<std::string, float> simulate(Particles& particles, AABB sim_area, float dt, vec2f gravity, int numPressureIters, int numParticleIters, float overRelaxation, bool compensateDrift);
    void draw(AABB area, sf::RenderTarget &window,
              std::unordered_map<eCellTypes, Color> color_table = {
              {eCellTypes::Air, Color(0, 0, 50)},
              {eCellTypes::Solid, Color(100, 100, 100)},
              {eCellTypes::Fluid, Color(70, 100, 220)}});
    Fluid(float cell_size, int width, int height);
};
#endif
