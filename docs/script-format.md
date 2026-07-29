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

### `accel_sq_color_thresh` — Acceleration Color Threshold

```
accel_sq_color_thresh  <value>
```

Threshold in acceleration-squared units used to determine the red channel color intensity of particles. The color mapping uses `cbrt(acc_sq / threshold)` to compute the red intensity, producing a blue-to-red gradient based on how strongly each particle is being accelerated relative to this threshold.

- Higher values shift more particles toward blue (cold/low acceleration appearance)
- Lower values shift more particles toward red (hot/high acceleration appearance)

**Default:** 1000000.0

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

## Gravity Method

### `Gravity` — Force Calculation Method

```
Gravity  <Octree | P2P>
```

Selects how gravitational forces are computed:

- **`Octree`** — Barnes-Hut tree algorithm. O(N log N) complexity. Required for large particle counts (>1000). Accuracy controlled by `BH_Opening_Theta`.
- **`P2P`** — Direct particle-to-particle (all-pairs). O(N^2) complexity. Exact forces. Only practical for small N (< ~1000 bodies).

**Default:** Octree

---

## Display

### `Display` — Render Resolution

```
Display  <width>  <height>
```

Sets the resolution of the OpenGL rendering viewport (and video recording output, if enabled). This controls the pixel dimensions of the rendered content, not including any window border or title bar added by the operating system.

**Examples:**
```
Display  1920  1080    # Full HD
Display  1280  720     # HD (default)
Display  3840  2160    # 4K
```

**Default:** 1280 720

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

When enabled (`1`), records the simulation output to an MP4 video file (`output.mp4`) using FFmpeg. Each rendered frame is captured and encoded at 30 FPS. Frames are only recorded while the simulation is running — pausing the simulation pauses recording.

**Default:** 0 (disabled)

---

### `End_Time` — Simulation End Time

```
End_Time  <value>
```

Sets a simulation time at which the program will automatically exit. The value is in simulation time units (same units as `dt`). Useful for batch rendering videos of a fixed duration.

If not specified, the simulation runs indefinitely until manually closed (Escape or window close).

**Examples:**
```
End_Time  5.0       # Exit when simulation time reaches 5.0
End_Time  100.0     # Run for 100 time units
```

**Default:** None (runs indefinitely)

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

Sets the initial 3D position of the camera in **graphical units** (i.e., after `DisplayScale` is applied). The camera looks toward the look-at point (default: origin `(0, 0, 0)`, configurable via `Camera_lookAt`).

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
- I/K — shift camera look-at point up/down (Y axis)

---

### `Camera_lookAt` — Camera Look-At Point

```
Camera_lookAt  <x>  <y>  <z>
```

Sets the point in space that the camera looks at and orbits around. The camera position (`Camera`) is interpreted as an offset relative to this point for computing spherical orbit angles. Can appear before or after `Camera` in the script.

**Examples:**
```
Camera_lookAt  0.0  0.0  0.0        # Look at the origin (default)
Camera_lookAt  1500.0  0.0  0.0     # Look at a point offset in +x (e.g., midpoint between two galaxies)
Camera_lookAt  0.0  100.0  0.0      # Look at a point above the orbital plane
```

**Default:** `0.0  0.0  0.0` (origin)

---

### `Camera_Orbit` — Automatic Camera Orbiting

```
Camera_Orbit  <theta_per_frame>
```

Enables automatic camera orbiting around the Y axis. Each simulation frame, the camera rotates by the specified theta angle (in radians). The orbit only advances when the simulation is running (not when paused).

| Parameter | Description |
|---|---|
| `theta_per_frame` | Angle in radians to rotate the camera per simulation frame. Positive values orbit counter-clockwise when viewed from above (+Y). |

**Examples:**
```
Camera_Orbit  0.01     # Moderate orbit speed (~6.3 seconds per revolution at 60 fps)
Camera_Orbit  0.005    # Slow orbit
Camera_Orbit  -0.01    # Orbit in the opposite direction
```

**Default:** Disabled (no automatic orbiting)

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

The dark matter halo applies a cored isothermal sphere potential: `a_halo = v_c^2 * r / (r^2 + r_c^2)`, centered on the mass-weighted barycenter of that system's particles (recomputed every derivative evaluation, not fixed to the central body). Every particle feels its own system's halo plus the halo of every other system, so multiple galaxies interact through their halos as well as their particles. The halo is a rigid analytic background — it never deforms, is never tidally stripped, and has no truncation radius; see `docs/dark-matter-halo.md` for what this does and does not model.

**Example (Milky-Way-like galaxy with dark matter halo, disc in x-z plane):**
```
GalaxyDisc  0   0.0 0.0 0.0   0.0 0.0 0.0   0.0 1.0 0.0   1.0e7 0.5 250.0 25.0 0.1  200.0 50.0
```

---

### `SphericalUniverse` — Procedural Spherical Distribution

```
SphericalUniverse  <system>  <posX> <posY> <posZ>  <velX> <velY> <velZ>  <total_mass> <radius> <H>  <halo_circular_velocity> <halo_core_radius>
```

Generates a uniform spherical distribution of particles with Hubble flow initial velocities. Each particle receives an outward radial velocity proportional to its distance from the center (v = H * r), mimicking cosmological expansion. Useful for simulating large-scale structure formation, collapsing gas clouds, or cosmological initial conditions.

| Parameter | Description |
|---|---|
| `system` | System index (0-based) |
| `posX`, `posY`, `posZ` | Position of the sphere center |
| `velX`, `velY`, `velZ` | Bulk velocity of the entire sphere |
| `total_mass` | Total mass distributed equally among all particles |
| `radius` | Radius of the spherical distribution. Particles are placed uniformly in volume (r scales as `rand^(1/3)`) with a small Gaussian perturbation (sigma=10 code units) |
| `H` | Hubble parameter. Each particle receives a radial outward velocity v = H * r, where r is its distance from the sphere center. This establishes Hubble flow (expansion proportional to distance). The critical value H_crit = sqrt(2*G*M/R^3) gives marginal unbinding; H < H_crit produces a bound (closed) system, H > H_crit produces an unbound (open) system. Set to 0 for no initial expansion (static sphere). |
| `halo_circular_velocity` | Dark matter halo circular velocity (same model as GalaxyDisc). Set to `0.0` for no halo |
| `halo_core_radius` | Dark matter halo core radius. Irrelevant if `halo_circular_velocity` is 0 |

**Example (100k-body expanding sphere with Hubble flow for structure formation):**
```
SphericalUniverse  0   0.0 0.0 0.0   0.0 0.0 0.0   5.0e7 200.0 2.83  0.0 1.0
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

Gravity         Octree

Display         1920  1080
DataLog         0
RecordVideo     1
Camera_Orbit    0.005

N_SystemBodies  30000  10000

Camera          0.0  500.0  500.0
Camera_lookAt   150.0  0.0  -250.0

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

---

## Command Line Usage

```
N-Body [script_file]
```

| Argument | Description |
|---|---|
| `script_file` | Path to a `.sim` script file. Default: `scripts/default.sim` |

**Examples:**
```bash
N-Body                                   # Run default script indefinitely
N-Body scripts/M51.sim                   # Run M51 script
```
