
![Demo](https://github.com/Epim3dium/flip-fluid-parallel/blob/fa88775d0ce31b8534089d728dee8634e7e5df3f/showoff.gif)
# Symulacja Cieczy Real-Time

## Streszczenie

Celem projektu jest porównanie trzech implementacji tego samego podejścia do symulacji cieczy w 2D, a konkretniej dwu fazowej symulacji wody i powietrza oraz tego jak wchodzą w interakcje. Aby pozwolić na eksperymentację z symulacją i obserwacje ciekawych zjawisk zachodzących w przekroju poprzecznym wody symulacja ma się odbywać w czasie rzeczywistym.


#### Podstawy teoretyczne
  * Algorytm symulacji FLIP (rozszerzenie PIC)
    * https://pub.ista.ac.at/group_wojtan/projects/2016_Ferstl_NBFLIP/nbflip.pdf
    * https://people.csail.mit.edu/kuiwu/gvdb_sim.html
    * https://en.wikipedia.org/wiki/Particle-in-cell
    * https://matthias-research.github.io/pages/tenMinutePhysics/18-flip.pdf
  * Podstawy matematyczne
    * https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations
    * https://pl.wikipedia.org/wiki/Ca%C5%82kowanie_numeryczne
    * https://pl.wikipedia.org/wiki/Metoda_Gaussa-Seidla
  * Detekcja kolizji (grid lookup):
    * https://www.gorillasun.de/blog/particle-system-optimization-grid-lookup-spatial-hashing/
#### Sprzęt
CPU: Intel i5-12600k \
GPU: RTX 4070 Super \
RAM: 32GB

## Opis problemu
Symulacja typu FLIP jest połączeniem eulerycznej symulacja płynów oraz podejściu symulacji cieczy za pomocą cząsteczek. Poprzez używanie dużej ilości małych cząstek zachowany jest niski współczynnik lepkości cieczy, oraz stabilność która jest dostarczona przez część euleryczną. Dla ułatwienia implementacji zakładam nieściśliwość wody. Aby rozwiązać system równań związanych z nieściśliwością i transfer prędkości z 'grid' do cząsteczek używam metody gaussa-siedela (wielu iteracji metody relaksacyjnej). Głównym wąskim gardłem takiej symulacji jest detekcja i rozstrzygnięcie kolizji między cząsteczkami więc na tym głównie się skupiam w moich implementacjach równoległych.
### Zarys symulacji
Cała pętla symulacji składa się z:
  - Symulacja cząsteczek
    * Akceleracja (grawitacja)
    * Scałkowanie prędkości i pozycji (dodanie akceleracji do prędkości i prędkości do pozycji)
    * Ograniczenie pozycji cząsteczek do pola symulacji
    * **Detekcja i rozwiązanie kolizji między cząsteczkami** (najbardziej czasochłonne)
  - Transfer prędkości z cząsteczek do macierzy (zapis kopii)
  - Obliczenie nowych prędkości które wymuszają nieściśliwość cieczy
  - Dodanie zmian w prędkości do cząsteczek (porównanie do kopii)
Z moich pomiarów wynika, iż **Detekcja i rozwiązanie kolizji między cząsteczkami** jest najbardziej czasochłonnym krokiem, więc na tym skupiłem swoje wysiłki w urównolegleniu.

### Grid-Lookup
Jedną z podstawowych podejść do optymalizacji detekcji kolizji wielu (okrągłych) cząsteczek o tych samych rozmiarach to grid-lookup. Polega na podzieleniu obszaru symulacji na siatkę gdzie każda komórka jest o rozmiarze $2R$x$2R$ gdzie $R$ to promień cząsteczek. Następnie porównywane są cząsteczki tylko z sąsiadujących komórek, w pionie, w poziomie i w ukosie, co optymistycznie znacznie przyspiesza wyszukiwanie kolizji. Szczególnie w obecnym przypadku, gdzie ciecz jest nieściśliwa rozłożenie cząsteczek powinno być równomierne, lecz w przypadku pesymistycznym złożoność tego algorytmu to nadal $ O(n^2) $

### Spis zaimplementowanych algorytmów (alogrytmów kolizji)

| Lp   | KOD   | Kategoria                    | Przeznaczenie                                                                          |
| ---- | ----- | ----------------------       | -----------------------------------------------------                                  |
| 1    | SEQ   | Algorytm Sekwencyjny         | Podstawowa implementacja algorytmu Grid-Lookup                                         |
| 2    | OMP   | Algorytm Równoległy (OpenMP) | Implementacja Grid-Lookup pozwalajaca na zrównoleglenie detekcji i rozwiązania kolizji |
| 3    | CUDA  | Algorytm Równoległy (CUDA)   | Przeniesienie wszystkich obliczeń związanych z cząsteczkami na kartę graficzną         |

## Implementacja
Ponieważ jest to bardzo duży projekt postaram się zwrzeć tylko najważniejsze fragmenty
Dla całego kodu proszę się skierować do [[https://github.com/Epim3dium/flip-fluid-parallel]]
```
fluid-sim [OPCJE]
Opcje:
  -w [szerokość] 
  -h [wysokość] 
  -c [rozmiar komórki FLIP] 
  -r [promień cząsteczki] 
  -n [ilość cząsteczek]
  -i [ilość iteracji solvera cząsteczek]
  -f [ilość iteracji solvera nieściślowości]
  -x [wsp. relaksacji]
  -d [0/1 czy rozpychać ściśnione cząsteczki]
  -t [czas raportowania w sekundach]
```
Wyjście w formacie json:
```
{
	"particle radius" : "5",
	"particle count" : "1000",
	"num of particle iters" : "32",
	"num of fluid iters" : "32",
	"overrelaxation" : "1.9",
	"raporting inter" : "0.125",
	"measurements" : [
		{
			"fluid::density" : "1.22811e-05",
			"fluid::incompressibility" : "6.32468e-05",
			"fluid::transfer1" : "4.48991e-05",
			"fluid::transfer2" : "4.35779e-05",
			"particles::accelerate" : "9.69806e-06",
			"particles::collide" : "0.001483",
			"particles::collide::assign" : "0.000765755",
			"particles::collide::cleanup" : "0.000117517",
			"particles::collide::compare" : "0.000307767",
			"particles::constraint" : "5.2476e-05",
			"particles::integrate" : "1.49574e-05",
			"FPS" : "343.032"
		}, ...]
}
```
#### Struktura cząsteczek
```
struct Particles {
    uint32_t max_particles_count;
    float radius;
    float diameter; // radius*2
    vec2f *position;
    vec2f *velocity;
    vec2f *acceleration;
};
```
### Głowna pętla symulacji
```
std::map<std::string, float> Fluid::simulate(Particles& particles, AABB sim_area, float dt, vec2f gravity, int numPressureIters, int numParticleIters, float overRelaxation, bool compensateDrift) {
    auto sdt = dt / numSubSteps;
    std::map<std::string, float> bench;

    sim_area.setSize(sim_area.size() - vec2f(m_cell_size, m_cell_size)*2.f);

    float col_time = 0.f;
    float fluid_time = 0.f;
    for (int step = 0; step < numSubSteps; step++) {
        Stopwatch local_stop;
        {
            ParticleSolveBlock solv(particles);
            for(int i = 0; i < numParticleIters; i++) {
                accelerate(particles, gravity);
                bench["particles::accelerate"] += local_stop.restart();
                integrate(particles, sdt / (float)numParticleIters);
                bench["particles::integrate"] += local_stop.restart();
                constraint(particles, sim_area);
                bench["particles::constraint"] += local_stop.restart();
                auto results = collide(particles, sim_area);
                for(auto [key, v] : results) bench[key] += v;
                bench["particles::collide"] += local_stop.restart();
            }
        }
        local_stop.restart();
        transferVelocities(true, 1.f, particles);
        bench["fluid::transfer1"] += local_stop.restart();
        updateParticleDensity(particles);
        bench["fluid::density"] += local_stop.restart();
        solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
        bench["fluid::incompressibility"] += local_stop.restart();
        transferVelocities(false, flipRatio, particles);
        bench["fluid::transfer2"] += local_stop.restart();
    }
    return bench;
}
```
### SEQ - Algorytm sekwencyjny (Grid lookup)
Podstawowa implementacja Grid Lookup
```
void processCollision(Particles& particles, int i, int ii) {
    auto diff = particles.position[ii] - particles.position[i];
    const float min_dist = particles.radius * 2;
    if(qlen(diff) < min_dist * min_dist && qlen(diff) > 1e-10f) {
        auto l = length(diff);
        auto n = normal(diff);
        static const float damping = 0.7f;
        auto c = (min_dist - l) * 0.5f * damping;
        particles.position[i] -= n * c;
        particles.position[ii] += n * c;
        resolveVelocities(particles, i, ii, -n);
    }
}

void compareWithNeighbours(Particles& particles, int col, int row, int max_segs_rows, const std::vector<CompactVec>& grid) {
    auto& comp_vec1 = grid[row*max_segs_rows + col];
    if(comp_vec1.size() == 0) return;
    #define offset_grid(dirx, diry) &grid[(row+dirx) * max_segs_rows + col+diry]
    std::array<const CompactVec*, 4U> comp_vecs = {offset_grid(1, 0), offset_grid(0, 1), offset_grid(1, 1),offset_grid(1,-1)};

    for(int i = 0; i < comp_vec1.size(); i++) {
        auto idx1 = comp_vec1[i];
        for(int ii = i + 1; ii < comp_vec1.size(); ii++) {
            auto idx2 = comp_vec1[ii];
            processCollision(particles, idx1, idx2);
        }
        for(auto other : comp_vecs) {
            for(auto ii = 0; ii < other->size(); ii++) {
                auto idx2 = (*other)[ii];
                processCollision(particles, idx1, idx2);
            }
        }
    }
}
std::map<std::string, float> collide(Particles& particles, AABB sim_area) {
    std::map<std::string, float> result;
    const uint32_t max_segs_cols = sim_area.size().x / particles.diameter + 1 + 2;
    const uint32_t max_segs_rows = sim_area.size().y / particles.diameter + 1 + 2;

    if(col_grid.size() != max_segs_rows * max_segs_cols) {
        col_grid = std::vector<CompactVec>(max_segs_rows*max_segs_cols);
    }
    auto max_dim = std::max(max_segs_cols, max_segs_rows);
    std::unordered_set<uint32_t> active_containers;
    Stopwatch stop;
    int counter = 0;
    for(int i = 0; i < max_particle_count; i++) {
        uint32_t col = (particles.position[i].x - sim_area.min.x) / particles.diameter;
        uint32_t row = (particles.position[i].y - sim_area.min.y) / particles.diameter;
        if(col+1 >= max_segs_cols || row+1 >= max_segs_rows) {
            std::cerr << "out of range particle\n";
            continue;
        }
        auto& comp_vec = col_grid[(row+1) * max_segs_rows + col+1];
        comp_vec.push_back(i);
        active_containers.insert((row+1) * max_dim + col+1);
    }
    result["particles::collide::assign"] += stop.restart();

    for(auto i : active_containers) {
        auto row = i / max_dim;
        auto col = i % max_dim;
        compareWithNeighbours(particles, col, row, max_segs_rows, col_grid);
    }
    result["particles::collide::compare"] += stop.restart();

    for(auto container : active_containers) {
        auto row = container / max_dim;
        auto col = container % max_dim;
        col_grid[row * max_segs_rows + col].clear();
    }
    result["particles::collide::cleanup"] += stop.restart();
    return result;
}
```
### OMP - Algorytm równoległy, OpenMP
Aby pozwolić na równoległe przetwarzanie kolizje są przetwarzanie w 9 fazach. Podział na dziewięc faz jest kluczowy, poniważ gwarantuje brak race-conditions (każda komórka może przetworzyć wszystkich swoich sąsiadów bez obaw iż ktoś inny też ich przetwarza). Fragment reprezentujący podział na 9 faz to:
```
int checkerboard = (row % 3)*3 + col%3;
active_containers[checkerboard].insert((row+1) * max_dim + col+1);
```
Cały kod:
```
std::map<std::string, float> collide(Particles& particles, AABB sim_area) {
    std::map<std::string, float> result;
    const uint32_t max_segs_cols = sim_area.size().x / particles.diameter + 1 + 2;
    const uint32_t max_segs_rows = sim_area.size().y / particles.diameter + 1 + 2;
    std::vector<CompactVec>& col_grid = particles.col_grid;

    Stopwatch stop;
    if(col_grid.size() != max_segs_rows * max_segs_cols) {
        col_grid = std::vector<CompactVec>(max_segs_rows*max_segs_cols);
    }
    auto max_dim = std::max(max_segs_cols, max_segs_rows);
    std::unordered_set<uint32_t> active_containers[9];
    int counter = 0;
    for(int i = 0; i < particles.max_particle_count; i++) {
        uint32_t col = (particles.position[i].x - sim_area.min.x) / particles.diameter;
        uint32_t row = (particles.position[i].y - sim_area.min.y) / particles.diameter;
        if(col+1 >= max_segs_cols || row+1 >= max_segs_rows) {
            std::cerr << "out of range particle\n";
            continue;
        }
        auto& comp_vec = col_grid[(row+1) * max_segs_rows + col+1];
        comp_vec.push_back(i);
        int checkerboard = (row % 3)*3 + col%3;
        active_containers[checkerboard].insert((row+1) * max_dim + col+1);
    }
    result["particles::collide::assign"] += stop.restart();

    for(int i = 0; i < 9; i++) {
        std::vector<uint32_t> active(active_containers[i].begin(),active_containers[i].end());
        #pragma omp parallel for
        for(int i = 0; i < active.size(); i++) {
            auto row = active[i] / max_dim;
            auto col = active[i] % max_dim;
            compareWithNeighbours(particles, col, row, max_segs_rows, col_grid);
        }
    }
    result["particles::collide::compare"] += stop.restart();

    for(int i = 0; i < 9; i++) {
        for(auto container : active_containers[i]) {
            auto row = container / max_dim;
            auto col = container % max_dim;
            col_grid[row * max_segs_rows + col].clear();
        }
    }
    result["particles::collide::cleanup"] += stop.restart();
    return result;
}
```
### CUDA - Algorytm równoległy, CUDA
Pierwszą próbą było zportowanie kodu OpenMP do CUDY, lecz przy wyższej wartości iteracji solvera cząsteczek przesył danych do i z karty graficznej trwał za długo (16 razy na sekundę to dużo!). W zwiąsku z tym zdecydowałem się na przeniesienie wszyskich kalkulacji na cząsteczkach na kartę graficzną. Nowa struktura particles:
```
struct Particles {
    uint32_t max_particles_count;
    float radius = 3.0f;
    float diameter = radius*2;
    vec2f *position;
    vec2f *velocity;
    vec2f *acceleration;

    vec2f *gpu_position;
    vec2f *gpu_velocity;
    vec2f *gpu_acceleration;
    vec2f *gpu_position_dt;
    vec2f *gpu_velocity_dt;
};
```
Dodatkową metodą przyspieszenia jest dwu fazowego algorytmu. Przy pomocy dodatkowych tablic (gpu_position_dt i gpu_velocity_dt) gdzie przechowywane będą tylko delty, zmiany po kolizji aktualizowane będą atomowo. Dopiero po rozwiązaniu wszystkich kolizji będzie druga faza dodawania delt to pozycji i prędkości. (umieściłem tutaj jedynie ciekawe fragmenty zmian)
```
__global__ void compareWithNeighboursKernel(
    const vec2f* position, const vec2f* velocity,
    vec2f* position_dt, vec2f* velocity_dt,
    float radius, uint32_t* active, uint32_t* active_size, int max_segs_cols, int max_segs_rows, const CompactVec* grid) {
    int max_size = max_segs_rows * max_segs_cols;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= *active_size)
        return;
    auto row = active[i] / max_segs_cols;
    auto col = active[i] % max_segs_cols;
    auto& comp_vec1 = grid[row*max_segs_cols + col];
    if(comp_vec1.size == 0) 
        return;
    #define offset_grid(dirx, diry) &grid[(row+dirx) * max_segs_cols + col+diry]
    const CompactVec* comp_vecs[4] = {offset_grid(1, 0), offset_grid(0, 1), offset_grid(1, 1),offset_grid(1,-1)};

    for(int i = 0; i < comp_vec1.size; i++) {
        auto idx1 = comp_vec1.data[i];
        for(int ii = i + 1; ii < min(comp_vec1.size, COMP_SIZE); ii++) {
            auto idx2 = comp_vec1.data[ii];
            processCollision(position, velocity, position_dt, velocity_dt, idx1, idx2, radius);
        }
        for(int j = 0; j < 4U; j++) {
            auto other = comp_vecs[j];
            for(auto ii = 0; ii < min(other->size, COMP_SIZE); ii++) {
                auto idx2 = (*other).data[ii];
                processCollision(position, velocity, position_dt, velocity_dt, idx1, idx2, radius);
            }
        }
    }
}
__global__ void assignParticlesToGridKernel(
    const vec2f* position, float diameter,
    CompactVec* col_grid,
    vec2f sim_area_min,
    int max_segs_cols,
    int max_segs_rows,
    uint32_t* active, // Flattened grid of flags (per cell)
    uint32_t* active_size
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= max_particle_count) return;

    vec2f pos = position[i];
    uint32_t col = (pos.x - sim_area_min.x) / diameter;
    uint32_t row = (pos.y - sim_area_min.y) / diameter;

    if (col + 1 >= max_segs_cols || row + 1 >= max_segs_rows) return;

    int grid_idx = (row + 1) * max_segs_cols + (col + 1);
    CompactVec& comp_vec = col_grid[grid_idx];

    // Atomically get index to write into comp_vec
    uint32_t insert_idx = atomicAdd(&comp_vec.size, 1);
    if (insert_idx <= COMP_SIZE) {
        comp_vec.data[insert_idx] = i;
    }    
    bool old = (insert_idx != 0);  //only the first one insterts
    if(old == 1) return;

    int max_size = max_segs_cols * max_segs_rows;
    auto idx = atomicAdd(active_size, 1);
    active[idx] = grid_idx;
}

std::map<std::string, float> collide(Particles& particles, AABB sim_area) {
    std::map<std::string, float> result;
    const uint32_t max_segs_cols = sim_area.size().x / particles.diameter + 1 + 2;
    const uint32_t max_segs_rows = sim_area.size().y / particles.diameter + 1 + 2;

    CUDA_CALL(cudaMemset(particles.gpu_position_dt, 0, sizeof(vec2f) * max_particle_count));
    CUDA_CALL(cudaMemset(particles.gpu_velocity_dt, 0, sizeof(vec2f) * max_particle_count));
    Stopwatch stop;
    //allocate new if size changed
    if(col_grid_size != max_segs_rows * max_segs_cols) {
        if(col_grid_size != 0) {
            cudaFree(gpu_col_grid);
            cudaFree(gpu_active_idxs);
        }else {
            CUDA_CALL(cudaMalloc(&gpu_active_size, sizeof(uint32_t)));
            CUDA_CALL(cudaMemset(gpu_active_size, 0, sizeof(uint32_t)));
        }
        col_grid_size = max_segs_rows * max_segs_cols;
        CUDA_CALL(cudaMalloc(&gpu_active_idxs, sizeof(uint32_t) * col_grid_size));
        CUDA_CALL(cudaMalloc(&gpu_col_grid, sizeof(CompactVec) * col_grid_size));
        CUDA_CALL(cudaMemset(gpu_col_grid, 0, sizeof(CompactVec) * col_grid_size));
    }
    int threadsPerBlock = 1024;
    int blocks = (max_particle_count + threadsPerBlock - 1) / threadsPerBlock;
    assignParticlesToGridKernel<<<blocks, threadsPerBlock>>>(
        particles.gpu_position, particles.diameter,
        gpu_col_grid, 
        sim_area.min,
        max_segs_cols, max_segs_rows,
        gpu_active_idxs, gpu_active_size
    );
    cudaDeviceSynchronize();
    result["particles::collide::assign"] += stop.restart();

    compareWithNeighboursKernel<<<(col_grid_size + threadsPerBlock) / threadsPerBlock, threadsPerBlock>>>(
        particles.gpu_position,
        particles.gpu_velocity,
        particles.gpu_position_dt,
        particles.gpu_velocity_dt,
        particles.radius,
        gpu_active_idxs,
        gpu_active_size,
        max_segs_cols,
        max_segs_rows,
        gpu_col_grid);
    cudaDeviceSynchronize();
    result["particles::collide::compare"] += stop.restart();
    blocks = (max_particle_count + threadsPerBlock - 1) / threadsPerBlock;
    addDeltasKernel<<<blocks, threadsPerBlock>>>(particles.gpu_position, particles.gpu_position_dt, max_particle_count);
    addDeltasKernel<<<blocks, threadsPerBlock>>>(particles.gpu_velocity, particles.gpu_velocity_dt, max_particle_count);
    cudaDeviceSynchronize();

    blocks = (col_grid_size + threadsPerBlock - 1) / threadsPerBlock;
    // clearGrid<<<blocks, threadsPerBlock>>>(gpu_col_grid, col_grid_size);
    CUDA_CALL(cudaMemset(gpu_active_size, 0, sizeof(uint32_t)));
    CUDA_CALL(cudaMemset(gpu_col_grid, 0, sizeof(CompactVec) * col_grid_size));
    CUDA_CALL(cudaMemset(gpu_active_idxs, 0, sizeof(uint32_t) * col_grid_size));
    result["particles::collide::cleanup"] += stop.restart();
    return result;
}

```
Oraz potrzebna dodatkowa struktura do przesyłu danych do i z karty graficznej na w bloku solvera cząsteczek. Potrzeba ParticleSolveBlock w głownej pętli:
```
ParticleSolveBlock::ParticleSolveBlock(Particles& p) : particles(p) {
    CUDA_CALL(cudaMemcpy(particles.gpu_velocity, particles.velocity, sizeof(vec2f) * max_particle_count, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(particles.gpu_position, particles.position, sizeof(vec2f) * max_particle_count, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(particles.gpu_acceleration, particles.acceleration, sizeof(vec2f) * max_particle_count, cudaMemcpyHostToDevice));
}
ParticleSolveBlock::~ParticleSolveBlock() {
    CUDA_CALL(cudaMemcpy(particles.velocity, particles.gpu_velocity, sizeof(vec2f) * max_particle_count, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(particles.position, particles.gpu_position, sizeof(vec2f) * max_particle_count, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(particles.acceleration, particles.gpu_acceleration, sizeof(vec2f) * max_particle_count, cudaMemcpyDeviceToHost));
}
```
## Testowanie wydajności
Ponieważ symulacja ma być real-time, jedynym wymogiem jest zachowanie 60 klatek na sekundę. Dla każdego z podejść przygotowałem serię pomiarów aby sprawdzić jak się zachowują pod innymi typami prac. \
W każdym z testów ciecz będzie startować w kolumnie a pomiary będą się odbywać do uzyskania statycznej tafli wody.
### Parametry testów
Dla każdego testu: 
  * $\text{wsp. relaksacji} = 1.9$
  * $\text{rozmiar komórki siatki cieczy} = \text {promień cząsteczki} * 2$
  * $\text{Liczba iteracji obliczeń cieczy} = 32$
(L.I. oznacza Liczbę iteracji)
| Nr. Testu   | Liczba cząsteczek   | Promień cząsteczki   | L.I. obliczeń cząsteczek   |
| ----------- | ------------------- | -------------------- | -------------------------- |
| 1           | 1000                | 5                    | 32                         |
| 2           | 2500                | 3                    | 32                         |
| 3           | 5000                | 2.5                  | 16                         |
| 4           | 5000                | 2.5                  | 32                         |
| 5           | 10000               | 2                    | 16                         |

W wynikach testów zawarłem czas poświęcony na obliczanie kolizji (im mniej tym lepiej), ponieważ to tam różnią się implementacje oraz miarę klatek na sekundę (im więcej tym lepiej) - kluczowe do oceny czy dany algorytm pozwala na interakcję w czasie rzeczywistym.

### Test 1
$\text{Liczba cząsteczek} = 1000$\
$\text{Promień cząsteczki} = 5$\
$\text{Liczba iteracji solvera kolizji} = 32$\
![11](https://github.com/Epim3dium/flip-fluid-parallel/blob/6a6c75a61c561038be48767e7c9564065c5b1ddc/media/11.png)
![12](https://github.com/Epim3dium/flip-fluid-parallel/blob/6a6c75a61c561038be48767e7c9564065c5b1ddc/media/12.png)
### Test 2
$\text{Liczba cząsteczek} = 2500$\
$\text{Promień cząsteczki} = 3$\
$\text{Liczba iteracji solvera kolizji} = 32$\
![21](https://github.com/Epim3dium/flip-fluid-parallel/blob/6a6c75a61c561038be48767e7c9564065c5b1ddc/media/21.png)
![22](https://github.com/Epim3dium/flip-fluid-parallel/blob/6a6c75a61c561038be48767e7c9564065c5b1ddc/media/22.png)
### Test 3
$\text{Liczba cząsteczek} = 5000$\
$\text{Promień cząsteczki} = 2.5$\
$\text{Liczba iteracji solvera kolizji} = 16$ \
![31](https://github.com/Epim3dium/flip-fluid-parallel/blob/6a6c75a61c561038be48767e7c9564065c5b1ddc/media/31.png)
![32](https://github.com/Epim3dium/flip-fluid-parallel/blob/6a6c75a61c561038be48767e7c9564065c5b1ddc/media/32.png)
### Test 4
$\text{Liczba cząsteczek} = 5000$\
$\text{Promień cząsteczki} = 2.5$\
$\text{Liczba iteracji solvera kolizji} = 32$\
![41](https://github.com/Epim3dium/flip-fluid-parallel/blob/6a6c75a61c561038be48767e7c9564065c5b1ddc/media/41.png)
![42](https://github.com/Epim3dium/flip-fluid-parallel/blob/6a6c75a61c561038be48767e7c9564065c5b1ddc/media/42.png)
### Test 5
$\text{Liczba cząsteczek} = 10000$\
$\text{Promień cząsteczki} = 2$\
$\text{Liczba iteracji solvera kolizji} = 16$\
![51](https://github.com/Epim3dium/flip-fluid-parallel/blob/6a6c75a61c561038be48767e7c9564065c5b1ddc/media/51.png)
![52](https://github.com/Epim3dium/flip-fluid-parallel/blob/6a6c75a61c561038be48767e7c9564065c5b1ddc/media/52.png)

## Wnioski
Łatwość zrównoleglenia kalkulacji kolizji sprawia iż zarówno implementacja OpenMP jak i CUDA kompletnie dominuje w pierwszych dwóch przypadkach testowych i dopiero wraz ze wzrostem ilości cząsteczek, różnica wydajności implementacji CUDA oraz OpenMP zaczyna być znacząco widoczna. \
Implementacja sekwencyjna przestaje być interaktywna (klatki na sekunde spadają poniżej 60) już w przypadku testu numer 2. Oznacza to iż nie jest to rozwiązanie chociażby dla umiarowanego rozmiaru symulacji.
W testach 4' oraz 5' implementacje OpenMP ledwo pozwala na interakcje. Nie jest więc dobrym wyborem dla bardziej wymagających symulacji. \
Porównując wyniki CUDA z testów 3 oraz 4 można stwierdzić iż nie różnią się one znacznie, pomimo iż w teście 4' wykonywane są dwa razy więcej iteracji. Można z tego wywnioskować, iż głównym wąskim gardłem dla implementacji CUDA jest czas przesyłu danych z procesora na kartę graficzną i z powrotem. Sprawia to, że jest to idealne podejście dla symulacji większych/bardziej dokładnych symulacji  cieczy (o większej wartości liczby iteracji solvera). \
Kolejnym etapem projektu byłoby zrównoleglenie części odpowiedzialnej za rozwiązywanie równania nieściśliwości. 
