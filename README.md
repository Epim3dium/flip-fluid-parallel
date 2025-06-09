# FLIP Fluid simulation
### Abstract

The aim of the project is to compare three implementations of the same approach to fluid simulation in 2D - a two-phase simulation of water and air and their interaction. The simulation is to run in real time to allow observation of phenomena in the fluid cross-section.

## Demo

[Demo](https://github.com/Epim3dium/flip-fluid-parallel/blob/fa88775d0ce31b8534089d728dee8634e7e5df3f/showoff.gif)

### Theoretical basis

- **Mathematical basis**.
  - Navier-Stokes equations](https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations)
  - Numerical integration](https://pl.wikipedia.org/wiki/Ca%C5%82kowanie_numeryczne)
  - Gauss-Seidel method](https://pl.wikipedia.org/wiki/Metoda_Gaussa-Seidla)
- **FLIP simulation algorithm (PIC extension)**.
  - Particle-in-cell](https://en.wikipedia.org/wiki/Particle-in-cell)
  - FLIP - ten minute physics](https://matthias-research.github.io/pages/tenMinutePhysics/18-flip.pdf)
- **Collision detection (grid lookup)**.
  - [Grid Lookup and Spatial Hashing](https://www.gorillasun.de/blog/particle-system-optimization-grid-lookup-spatial-hashing/)

### Hardware

- **CPU**: Intel i5-12600k  
- **GPU**: RTX 4070 Super  
- **RAM**: 32GB  

### Problem description

FLIP simulation combines euleric and particle-based approaches. The use of many small particles keeps the simulation low viscosity and stable. The incompressibility of water is assumed. A key aspect is **detection and collision resolution**, a major performance bottleneck. Implementations focus on optimizing this step.

### Simulation outline

1 **Particle simulation**.
   - Gravity
   - Integrating velocity and position
   - Position restriction to the simulation field
   - **Collision detection and resolution**.
2 **Velocity transfer**.
3 **Solving incompressibility**.
4 **Particle velocity update**.

### Grid-Lookup

The simulation field is divided into a grid of `2R x 2R` cells, where `R` is the particle radius. Particles are compared only within neighboring cells, which optimizes complexity. In the pessimistic case, it still amounts to `O(n^2)`.

### List of implemented algorithms (collision detection).

| Lp   | Code   | Category                    | Purpose                                                                     | Notes   |
| ---- | ------ | --------------------------- | --------------------------------------------------------------------------- | ------- |
| 1    | SEQ    | Sequential Algorithm        | Basic Grid-Lookup Implementation                                            |
| 2    | OMP    | Parallel (OpenMP)           | Grid-Lookup with Parallel Collision Detection and Resolution                |
| 3    | CUDA   | Parallel (CUDA)             | All calculations moved to GPU                                               | 
