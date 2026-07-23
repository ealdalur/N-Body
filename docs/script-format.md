# Simulation Script Format (.sim)

The N-Body simulator uses plain-text script files (`.sim`) to define all simulation parameters and initial conditions. Scripts replace all hard-coded constants, allowing different simulations to be configured without recompiling.

## General Syntax

- One command per line
- Comments start with `#` and extend to the end of the line
- Blank lines are ignored
- Whitespace separates the command keyword from its arguments
- Commands are case-sensitive
- Order matters: `N_SystemBodies` must appear before any initial condition commands (`Body`, `GalaxyDisc`, `SphericalUniverse`)

---

## Simulation Parameters

### `G` — Gravitational Constant

```
G  <value>
```

Sets the gravitational constant used in force calculations. The value depends on your unit system:

| Unit system | G value |
|---|---|
| Galactic (60 pc, 10^4 M_sun, 1 km/s) | 1.0 |
| Solar system (AU, solar masses, years) | 39.4784176044 (= 4pi^2) |

**Default:** 1.0

---

### `FDE` — Force of Dark Energy

```
FDE  <value>
```

Coefficient for a cosmological repulsive force that scales linearly with distance from the origin. When non-zero, each particle receives an outward radial acceleration proportional to `FDE * position`. Useful for confining spherical distributions or simulating expanding-universe effects.

- `0.0` = no dark energy (gravity only)
- `1.0` = moderate confining force

**Default:** 0.0

---

### `dt` — Time Step

```
dt  <value>
```

The fixed integration time step. Smaller values increase accuracy but slow the simulation. The required value depends on the fastest orbital period in your system:

| Simulation type | Typical dt |
|---|---|
| Galaxy (G=1 units) | 0.0005 |
| Solar system (AU/yr units) | 0.0001 |

**Rule of thumb:** Aim for at least 500-1000 time steps per shortest orbital period.

**Default:** 0.0005

---

### `r_soft` — Gravitational Softening Radius

```
r_soft  <value>
```

Softening length added to gravitational force calculations to prevent divergence at close encounters. The pairwise force denominator becomes `(r^2 + r_soft^2)^(3/2)` instead of `r^3`.

- Large values (e.g., 0.1) smooth forces at short range — appropriate for galaxy simulations where particles represent many stars.
- Small values (e.g., 1e-6) preserve point-mass behavior — appropriate for solar system simulations.

**Default:** 0.1

---

### `BH_Opening_Theta` — Barnes-Hut Opening Angle

```
BH_Opening_Theta  <value>
```

The opening angle parameter (theta) for the Barnes-Hut octree gravity algorithm. Controls the trade-off between accuracy and speed:

- `0.0` = exact N-body (never approximates, equivalent to P2P)
- `0.5` = balanced accuracy and performance (recommended)
- `1.0` = aggressive approximation (faster, less accurate)

Force error scales approximately as theta^2. Only relevant when `Gravity` is set to `Octree`.

**Default:** 0.5

---

### `DisplayScale` — Rendering Scale Factor

```
DisplayScale  <value>
```

Multiplies particle positions at render time only. Physics is unaffected — this purely controls how spread out particles appear on screen. Useful when simulation distances are much smaller or larger than the camera controls expect (which are tuned for separations of 100-500 graphical units).

- `1.0` = 1 simulation unit = 1 graphical unit
- `4.0` = 1 simulation unit = 4 graphical units (e.g., solar system in AU)
- `10.0` = 1 simulation unit = 10 graphical units

The camera position (set via `Camera`) is specified in graphical units, so adjust it accordingly when changing this value.

**Default:** 1.0

---

## Solver and Method

### `Solver` — Integration Method

```
Solver  <LeapFrog | RK4>
```

Selects the numerical integration scheme:

- **`LeapFrog`** — Velocity Stormer-Verlet (symplectic). Conserves energy over long timescales. Recommended for most simulations, especially those running for many orbital periods.
- **`RK4`** — Classical 4th-order Runge-Kutta. Higher instantaneous accuracy per step but does not conserve energy — orbits will slowly spiral inward over long runs.

**Default:** LeapFrog

---

### `Gravity` — Force Calculation Method

```
Gravity  <Octree | P2P>
```

Selects how gravitational forces are computed:

- **`Octree`** — Barnes-Hut tree algorithm. O(N log N) complexity. Required for large particle counts (>1000). Accuracy controlled by `BH_Opening_Theta`.
- **`P2P`** — Direct particle-to-particle (all-pairs). O(N^2) complexity. Exact forces. Only practical for small N (< ~1000 bodies).

**Default:** Octree

---

## Output Options

### `DataLog` — Binary Data Logging

```
DataLog  <0 | 1>
```

When enabled (`1`), writes binary state data (position magnitudes, velocity magnitudes, acceleration magnitudes) to a log file each time step. Useful for post-processing analysis.

**Default:** 0 (disabled)

---

### `RecordVideo` — Video Recording

```
RecordVideo  <0 | 1>
```

When enabled (`1`), records the simulation output to an MP4 video file (`output.mp4`) using FFmpeg. Each rendered frame is captured and encoded at 30 FPS.

**Default:** 0 (disabled)

---

## System Configuration

### `N_SystemBodies` — Number of Bodies Per System

```
N_SystemBodies  <count1>  [count2]  [count3]  ...
```

Defines how many particles belong to each gravitational system. Each space-separated integer defines one system. The total number of bodies is the sum of all counts.

A "system" is a group of particles that share a common central body (the first particle in the group) and optionally a dark matter halo. Multiple systems can interact gravitationally.

**Examples:**
```
N_SystemBodies  30000 10000    # Two systems: 30k + 10k = 40k total bodies
N_SystemBodies  40000          # One system with 40k bodies
N_SystemBodies  10             # One system with 10 bodies (e.g., solar system)
```

This command must appear before any initial condition commands.

---

## Camera

### `Camera` — Initial Camera Position

```
Camera  <x>  <y>  <z>
```

Sets the initial 3D position of the camera in **graphical units** (i.e., after `DisplayScale` is applied). The camera always looks toward the origin `(0, 0, 0)`.

The coordinate system for rendering:
- **x** — horizontal (positive = right on screen)
- **y** — vertical (positive = up from the orbital plane)
- **z** — depth (positive = toward bottom of screen when viewed from +y)

For a top-down view of the x-z plane, place the camera along the +y axis:
```
Camera  0.0  500.0  0.0      # 500 units above, looking down
Camera  0.0  500.0  500.0    # Angled view (45 degrees from above)
```

**Controls at runtime:**
- W/S — orbit camera up/down (phi angle)
- A/D — orbit camera left/right (theta angle)
- J/L — zoom in/out (radial distance)

---

## Initial Conditions

Initial condition commands define the starting state of all particles. They must appear after `N_SystemBodies`. Three types are available:

---

### `Body` — Single Explicit Body

```
Body  <system>  <posX> <posY> <posZ>  <velX> <velY> <velZ>  <mass>
```

Places a single particle with an explicitly specified state vector.

| Parameter | Description |
|---|---|
| `system` | System index (0-based) this body belongs to |
| `posX`, `posY`, `posZ` | Initial position in simulation units |
| `velX`, `velY`, `velZ` | Initial velocity in simulation units per time unit |
| `mass` | Mass of the body in simulation mass units |

Bodies are assigned to slots in order within their system. The first `Body` command for system 0 fills slot 0 (the central body), the second fills slot 1, and so on.

**Example (Sun and Earth in AU/solar-mass/year units):**
```
Body  0   0.0 0.0 0.0   0.0 0.0 0.0   1.0          # Sun
Body  0   1.0 0.0 0.0   0.0 0.0 -6.28  3.0e-6      # Earth
```

---

### `GalaxyDisc` — Procedural Galaxy Disc

```
GalaxyDisc  <system>  <posX> <posY> <posZ>  <velX> <velY> <velZ>  <normalX> <normalY> <normalZ>  <total_mass> <mass_fraction> <outer_radius> <inner_radius> <velocity_tolerance>  <halo_circular_velocity> <halo_core_radius>
```

Generates a flattened disc of particles with approximately circular orbits, representing a spiral galaxy. The disc plane is defined by the normal vector; particles orbit counter-clockwise when viewed from the direction the normal points.

| Parameter | Description |
|---|---|
| `system` | System index (0-based) |
| `posX`, `posY`, `posZ` | Position of the galaxy center |
| `velX`, `velY`, `velZ` | Bulk velocity of the entire galaxy |
| `normalX`, `normalY`, `normalZ` | Disc normal vector (defines orientation of the disc plane). Does not need to be unit length. Particles orbit counter-clockwise when viewed from the direction this vector points. Use `(0, 1, 0)` for a disc in the x-z plane |
| `total_mass` | Total mass of the central body (the galaxy core) |
| `mass_fraction` | Fraction of the total mass distributed among disc particles. `0.5` means disc particles collectively have half the central body's mass |
| `outer_radius` | Maximum radius of the disc |
| `inner_radius` | Minimum radius of the disc (creates a central hole) |
| `velocity_tolerance` | Random perturbation applied to orbital velocities. `0.0` = perfectly circular orbits, `0.1` = 10% random deviation |
| `halo_circular_velocity` | Asymptotic circular velocity of the dark matter halo. Controls how flat the rotation curve is at large radii. Set to `0.0` for no halo |
| `halo_core_radius` | Core radius of the dark matter halo (cored isothermal sphere). The halo density flattens inside this radius. Irrelevant if `halo_circular_velocity` is 0 |

The dark matter halo applies a cored isothermal sphere potential: `a_halo = v_c^2 * r / (r^2 + r_c^2)` centered on the central body.

**Example (Milky-Way-like galaxy with dark matter halo, disc in x-z plane):**
```
GalaxyDisc  0   0.0 0.0 0.0   0.0 0.0 0.0   0.0 1.0 0.0   1.0e7 0.5 250.0 25.0 0.1  200.0 50.0
```

---

### `SphericalUniverse` — Procedural Spherical Distribution

```
SphericalUniverse  <system>  <posX> <posY> <posZ>  <velX> <velY> <velZ>  <total_mass> <radius> <max_velocity>  <halo_circular_velocity> <halo_core_radius>
```

Generates a uniform spherical distribution of particles with random velocities. Useful for simulating collapsing gas clouds, globular clusters, or cosmological initial conditions.

| Parameter | Description |
|---|---|
| `system` | System index (0-based) |
| `posX`, `posY`, `posZ` | Position of the sphere center |
| `velX`, `velY`, `velZ` | Bulk velocity of the entire sphere |
| `total_mass` | Total mass distributed equally among all particles |
| `radius` | Radius of the spherical distribution. Particles are placed uniformly in volume (r scales as `rand^(1/3)`) with a small Gaussian perturbation |
| `max_velocity` | Maximum random velocity magnitude assigned to each particle. Directions are random (isotropic). Controls initial kinetic energy / "temperature" |
| `halo_circular_velocity` | Dark matter halo circular velocity (same model as GalaxyDisc). Set to `0.0` for no halo |
| `halo_core_radius` | Dark matter halo core radius. Irrelevant if `halo_circular_velocity` is 0 |

**Example (40k-body confined sphere):**
```
SphericalUniverse  0   0.0 0.0 0.0   0.0 0.0 0.0   2.1e7 210.0 200.0  0.0 1.0
```

---

## Complete Example

```
# Galaxy merger simulation
G               1.0
FDE             0.0
dt              0.0005
r_soft          0.1
BH_Opening_Theta  0.5

Solver          LeapFrog
Gravity         Octree

DataLog         0
RecordVideo     0

N_SystemBodies  30000  10000

Camera          0.0  500.0  500.0

GalaxyDisc  0   0.0 0.0 0.0       0.0 0.0 0.0       0.0 1.0 0.0    1.0e7 0.5 250.0 25.0 0.1  200.0 50.0
GalaxyDisc  1   300.0 0.0 -500.0  -150.0 0.0 0.0    0.0 1.0 0.0    6.0e6 0.5 125.0 25.0 0.1  155.0 25.0
```

---

## File Path Resolution

When the simulation loads a script, it searches for the file in three locations (in order):

1. The path as given (e.g., `scripts/default.sim`)
2. Prefixed with `../` (e.g., `../scripts/default.sim`)
3. Prefixed with `../../` (e.g., `../../scripts/default.sim`)

This allows the executable to find scripts regardless of whether it is run from the project root, the build directory, or a nested build configuration directory.
