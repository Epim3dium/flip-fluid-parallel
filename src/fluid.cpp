#include "fluid.hpp"
Fluid::Fluid(float cell_size, int width, int height) : m_cell_size(cell_size), m_width(width), m_height(height), m_num_cells(width * height){
    cell_color = std::vector<Color>(m_num_cells);
    cell_type = std::vector<eCellTypes>(m_num_cells, eCellTypes::Air);
    solid = std::vector<float>(m_num_cells, 1.f);
    particle_density = std::vector<float>(m_num_cells, 0.f);
    pressure = std::vector<float>(m_num_cells, 0.f);

    velocities = std::vector<vec2f>(m_num_cells, vec2f(0, 0));
    prev_velocities = std::vector<vec2f>(m_num_cells, vec2f(0, 0));
    velocities_diff = std::vector<vec2f>(m_num_cells, vec2f(0, 0));

    for(int i = 0; i < m_height; i++) {
        for(int j = 0; j < m_width; j++) {
            if(j == 0 || i == 0 || i == m_height - 1 || j == m_width - 1)
                solid[i + j * m_height] = 0.f;
        }
    }
}
void Fluid::updateParticleDensity(Particles& particles)
{
    int n = m_height;
    float h = m_cell_size;
    float h1 = 1.f / m_cell_size;
    float h2 = 0.5 * h;

    std::fill(particle_density.begin(), particle_density.end(), 0.f);

    for (int i = 0; i < max_particle_count; i++) {
        auto x = particles.position[i].x;
        auto y = particles.position[i].y;

        // x = std::clamp(x, h, (m_width - 1) * h);
        // y = std::clamp(y, h, (m_height - 1) * h);

        auto x0 = floorf((x - h2) * h1);
        auto tx = ((x - h2) - x0 * h) * h1;
        auto x1 = fmin(x0 + 1, m_width-2);
        
        auto y0 = floorf((y-h2)*h1);
        auto ty = ((y - h2) - y0*h) * h1;
        auto y1 = fmin(y0 + 1, m_height-2);

        auto sx = 1.0 - tx;
        auto sy = 1.0 - ty;

        if (x0 < m_width && y0 < m_height) particle_density[x0 * n + y0] += sx * sy;
        if (x1 < m_width && y0 < m_height) particle_density[x1 * n + y0] += tx * sy;
        if (x1 < m_width && y1 < m_height) particle_density[x1 * n + y1] += tx * ty;
        if (x0 < m_width && y1 < m_height) particle_density[x0 * n + y1] += sx * ty;
    }

    if (particleRestDensity == 0.0) {
        auto sum = 0.0;
        auto numFluidCells = 0;

        for (int i = 0; i < m_num_cells; i++) {
            if (cell_type[i] == eCellTypes::Fluid) {
                sum += particle_density[i];
                numFluidCells++;
            }
        }

        if (numFluidCells > 0)
            particleRestDensity = sum / numFluidCells;
    }
}
void Fluid::transferVelocities(bool toGrid, float flipRatio, Particles& particles)
{
    auto n = m_height;
    auto h = m_cell_size;
    auto h1 = 1.f / m_cell_size;
    auto h2 = 0.5 * h;

    if (toGrid) {

        prev_velocities = velocities;

        std::fill(velocities_diff.begin(), velocities_diff.end(), vec2f(0, 0));
        std::fill(velocities.begin(), velocities.end(), vec2f(0, 0));

        for (auto i = 0; i < m_num_cells; i++) 
            cell_type[i] = (solid[i] == 0.0 ? eCellTypes::Solid : eCellTypes::Air);

        for (auto i = 0; i < max_particle_count; i++) {
            auto x = particles.position[i].x;
            auto y = particles.position[i].y;
            auto xi = /* std::clamp( */floorf(x * h1)/* , 0.f, m_width - 1.f) */;
            auto yi = /* std::clamp( */floorf(y * h1)/* , 0.f, m_height - 1.f) */;
            auto cellNr = xi * n + yi;
            if (cell_type[cellNr] == eCellTypes::Air)
                cell_type[cellNr] = eCellTypes::Fluid;
        }
    }

    for (auto component = 0; component < 2; component++) {

        auto dx = component == 0 ? 0.0 : h2;
        auto dy = component == 0 ? h2 : 0.0;

        auto f = [&](int idx) -> float&{
            if(component == 0) {
                return velocities[idx].x;
            }
            return velocities[idx].y;
        };
        auto prevF = [&](int idx) -> float&{
            if(component == 0) {
                return prev_velocities[idx].x;
            }
            return prev_velocities[idx].y;
        };
        auto d= [&](int idx) -> float&{
            if(component == 0) {
                return velocities_diff[idx].x;
            }
            return velocities_diff[idx].y;
        };

        for (auto i = 0; i < max_particle_count; i++) {
            auto x = particles.position[i].x;
            auto y = particles.position[i].y;

            // x = std::clamp(x, h, (m_width - 1) * h);
            // y = std::clamp(y, h, (m_height - 1) * h);

            auto x0 = fmin(floorf((x - dx) * h1), m_width - 2);
            auto tx = ((x - dx) - x0 * h) * h1;
            auto x1 = fmin(x0 + 1, m_width-2);
            
            auto y0 = fmin(floorf((y-dy)*h1), m_height-2);
            auto ty = ((y - dy) - y0*h) * h1;
            auto y1 = fmin(y0 + 1, m_height-2);

            auto sx = 1.0 - tx;
            auto sy = 1.0 - ty;

            auto d0 = sx*sy;
            auto d1 = tx*sy;
            auto d2 = tx*ty;
            auto d3 = sx*ty;

            auto nr0 = x0*n + y0;
            auto nr1 = x1*n + y0;
            auto nr2 = x1*n + y1;
            auto nr3 = x0*n + y1;

            if (toGrid) {
                auto pv = particles.velocity[i].x;
                if(component == 1)
                    pv = particles.velocity[i].y;
                f(nr0) += pv * d0;  d(nr0) += d0;
                f(nr1) += pv * d1;  d(nr1) += d1;
                f(nr2) += pv * d2;  d(nr2) += d2;
                f(nr3) += pv * d3;  d(nr3) += d3;
            }
            else {
                auto offset = (component == 0 ? n : 1);
                auto valid0 = cell_type[nr0] != eCellTypes::Air || cell_type[nr0 - offset] != eCellTypes::Air ? 1.0 : 0.0;
                auto valid1 = cell_type[nr1] != eCellTypes::Air || cell_type[nr1 - offset] != eCellTypes::Air ? 1.0 : 0.0;
                auto valid2 = cell_type[nr2] != eCellTypes::Air || cell_type[nr2 - offset] != eCellTypes::Air ? 1.0 : 0.0;
                auto valid3 = cell_type[nr3] != eCellTypes::Air || cell_type[nr3 - offset] != eCellTypes::Air ? 1.0 : 0.0;

                auto v = particles.velocity[i].x;
                if(component == 1) 
                    v = particles.velocity[i].y;
                auto d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                if (d > 0.0) {

                    auto picV = (valid0 * d0 * f(nr0) + valid1 * d1 * f(nr1) + valid2 * d2 * f(nr2) + valid3 * d3 * f(nr3)) / d;
                    auto corr = (valid0 * d0 * (f(nr0) - prevF(nr0)) + valid1 * d1 * (f(nr1) - prevF(nr1))
                        + valid2 * d2 * (f(nr2) - prevF(nr2)) + valid3 * d3 * (f(nr3) - prevF(nr3))) / d;
                    auto flipV = v + corr;

                    if(component == 0) {
                        particles.velocity[i].x = (1.0 - flipRatio) * picV + flipRatio * flipV;
                    }else {
                        particles.velocity[i].y = (1.0 - flipRatio) * picV + flipRatio * flipV;
                    }
                }
            }
        }

        if (toGrid) {
            for (auto i = 0; i < velocities.size(); i++) {
                if (d(i) > 0.0)
                    f(i) /= d(i);
            }

            // restore solid cells

            for (auto i = 0; i < m_width; i++) {
                for (auto j = 0; j < m_height; j++) {
                    auto solid = (cell_type[i * n + j] == eCellTypes::Solid);
                    if (solid || (i > 0 && cell_type[(i - 1) * n + j] == eCellTypes::Solid))
                        velocities[i * n + j].x = prev_velocities[i * n + j].x;
                    if (solid || (j > 0 && cell_type[i * n + j - 1] == eCellTypes::Solid))
                        velocities[i * n + j].y = prev_velocities[i * n + j].y;
                }
            }
        }
    }
}
void Fluid::solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift) {
    std::fill(pressure.begin(), pressure.end(), 0.f);
    prev_velocities = velocities;

    auto n = m_height;
    auto cp = density * m_cell_size / dt;

    // for (auto i = 0; i < m_num_cells; i++) {
    //     auto u = u[i];
    //     auto v = v[i];
    // }

    for (auto iter = 0; iter < numIters; iter++) {

        for (auto i = 1; i < m_width-1; i++) {
            for (auto j = 1; j < m_height-1; j++) {

                if (cell_type[i*n + j] != eCellTypes::Fluid)
                    continue;

                auto center = i * n + j;
                auto left = (i - 1) * n + j;
                auto right = (i + 1) * n + j;
                auto bottom = i * n + j - 1;
                auto top = i * n + j + 1;

                auto sc = solid[center];
                auto sx0 = solid[left];
                auto sx1 = solid[right];
                auto sy0 = solid[bottom];
                auto sy1 = solid[top];
                auto s = sx0 + sx1 + sy0 + sy1;
                if (s == 0.0)
                    continue;

                auto div = velocities[right].x - velocities[center].x + 
                    velocities[top].y - velocities[center].y;

                if (particleRestDensity > 0.0 && compensateDrift) {
                    auto k = 1.0;
                    auto compression = particle_density[i*n + j] - particleRestDensity;
                    if (compression > 0.0)
                        div = div - k * compression;
                }

                auto dp = -div / s;
                dp *= overRelaxation;
                pressure[center] += cp * dp;

                velocities[center].x -= sx0 * dp;
                velocities[right].x += sx1 * dp;
                velocities[center].y -= sy0 * dp;
                velocities[top].y += sy1 * dp;
            }
        }
    }
}
void Fluid::simulate(Particles& particles, AABB sim_area, float dt, vec2f gravity, int numPressureIters, int numParticleIters, float overRelaxation, bool compensateDrift) {
    auto numSubSteps = 1;
    auto sdt = dt / numSubSteps;

    sim_area.setSize(sim_area.size() - vec2f(m_cell_size, m_cell_size)*2.f);

    for (int step = 0; step < numSubSteps; step++) {
        for(int i = 0; i < numParticleIters; i++) {
            accelerate(particles, gravity);
            integrate(particles, sdt / (float)numParticleIters);
            constraint(particles, sim_area);
            collide(particles, sim_area);
        }
        transferVelocities(true, 1.f, particles);
        updateParticleDensity(particles);
        solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
        transferVelocities(false, flipRatio, particles);
    }
}
void Fluid::draw(AABB area, sf::RenderTarget &window,
          std::unordered_map<eCellTypes, Color> color_table) {
    sf::Image img;
    img.create(m_width, m_height);
    for (int i = 0; i < m_height; i++) {
        for (int j = 0; j < m_width; j++) {
            auto type = cell_type[i + j * m_height];
            if (!color_table.contains(type)) {
                img.setPixel(j, m_height - i - 1, sf::Color(255, 0, 255));
            } else if (type == eCellTypes::Fluid) {
                auto point_density = particle_density[i + j * m_height];
                auto color = color_table.at(type);
                if(point_density > particleRestDensity) {
                    auto d = particleRestDensity / point_density  * 0.5 + 0.5f;
                    color.r *= d;
                    color.g *= d;
                    color.b *= d;
                }
                img.setPixel(j, m_height - i - 1, color);
            }else {
                img.setPixel(j, m_height - i - 1, color_table.at(type));
            }
        }
    }
    sf::Texture tex;
    tex.loadFromImage(img);
    sf::Sprite spr;
    spr.setTexture(tex);
    vec2f scale = vec2f(area.size().x / m_width, area.size().y / m_height);
    spr.setScale(scale);
    spr.setPosition(area.bl());
    window.draw(spr);
}
