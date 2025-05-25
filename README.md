# 3D Fluid Simulation with SPH & Spatial Hashing
*A real-time, interactive fluid simulator built with Three.js and WebGL, using Smoothed Particle Hydrodynamics (SPH) optimized via a 3D spatial grid.*


# Features

Physics Engine: Implements SPH for fluid dynamics (pressure, viscosity, gravity) with kernel functions (Poly6, Spiky, Viscosity).

Optimized Performance: Uses a 3D spatial grid to reduce neighbor-search complexity from O(n²) to O(n).

Interactive Controls:
- Adjust parameters (viscosity, gravity, particle count) via GUI.
- Use mouse hovering to interact with the fluid directly.
- Real-Time Rendering: WebGL-powered visualization with Three.js, including dynamic lighting and camera controls.
# Tech Stack

Frontend: JavaScript, Three.js, WebGL
Physics: SPH algorithms, spatial hashing, cubic spline kernels
# How It Works

Particles = Fluid: The system models fluid as thousands of interacting particles.
Grid Acceleration: Particles are sorted into a 3D grid to efficiently find neighbors.
Forces Calculated Per-Frame: Density, pressure, and viscosity forces update particle positions.

# Why It’s Cool

Educational: Clean code structure to learn SPH or spatial optimization.

Visual: Mesmerizing fluid behavior with tweakable parameters.

Performant: Handles 10k+ particles smoothly thanks to grid optimization.

For kernels and forces formula https://matthias-research.github.io/pages/publications/sca03.pdf
