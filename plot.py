import json
import matplotlib.pyplot as plt
import sys
import os


def displayFileData(filepath):
    # Read data from a JSON file
    with open(filepath, "r") as file:
        data = json.load(file)
    filename = os.path.basename(sys.argv[1])

    print("Simulation Information:")
    print(f"Particle Radius       : {data['particle radius']}")
    print(f"Particle Count        : {data['particle count']}")
    print(f"Particle Iterations   : {data['num of particle iters']}")
    print(f"Fluid Iterations      : {data['num of fluid iters']}")
    print(f"Overrelaxation Factor : {data['overrelaxation']}")
    print()

    iterations = list(range(len(data["measurements"])))
    fps = [float(m["FPS"]) for m in data["measurements"]]
    particle_time = [float(m["particle time"]) * 1000 for m in data["measurements"]]
    fluid_time = [float(m["fluid time"]) * 1000 for m in data["measurements"]]

    plt.figure(figsize=(10, 6))
    plt.plot(iterations, fps, marker='o', label='FPS', color='blue')
    plt.xlabel("Iteration")
    plt.ylabel("Frames per second")
    plt.title(f"FPS - {filename}")
    plt.ylim(0, 60)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    plt.plot(iterations, particle_time, marker='s', label='Particle Time [ms]', color='green')
    plt.plot(iterations, fluid_time, marker='^', label='Fluid Time [ms]', color='red')
    plt.xlabel("Iteration")
    plt.ylabel("Milliseconds")
    plt.title(f"Particle Time, and Fluid Time per Iteration - {filename}")
    plt.ylim(0, 60)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


displayFileData(sys.argv[1])
