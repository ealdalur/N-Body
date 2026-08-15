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

### `OrbitDiagnostic` — Two-Galaxy Orbit Logging

```
OrbitDiagnostic  <N>
```

For a run with two or more systems, appends one CSV row every `N` steps to `orbit_diagnostic.csv` (in the simulation's working directory) recording the two galaxies' barycentre and halo-centre positions/velocities, their separation, and the radial/tangential split of the relative velocity. Analyse with `scripts/analyze_orbit_diagnostic.py`, which reconstructs the specific orbital energy and compares the live orbit against the conservative analytic orbit — useful for detecting spurious orbital decay.

`N = 0` (or omitting the command) disables it entirely; there is no output and no per-step cost.

**Default:** 0 (disabled)

---

### `Info_Display` — On-Screen Information

```
Info_Display  <0 | 1>
```

When enabled (`1`), displays an information overlay on the rendered output showing FPS, simulation time (T), kinetic energy (KE), potential energy (PE), and total energy (E). When disabled (`0`), the overlay is hidden (useful for clean video recording).

**Default:** 1 (enabled)

---

### `InitializationTime` — Warmup / Settling Phase

```
InitializationTime  <value>
```

When greater than zero, the simulation starts at `t = -InitializationTime` with
the systems **dynamically isolated**, then transitions to normal coupled evolution
at `t = 0`. A value of `0` (the default) disables the feature entirely: the run
starts at `t = 0` and behaves exactly as it did before this was added.

Procedurally generated discs are placed in Jeans/Toomre equilibrium but still
carry particle-noise transients that produce odd, unphysical-looking structure for
the first fraction of a dynamical time. The warmup lets each system settle before
the interaction begins, so the encounter acts on a clean disc.

**During warmup (`t < 0`):**

| | Behaviour |
|---|---|
| Gravity | Only *within* each system. No cross-system particle forces. |
| Halo | Each system feels only its **own** halo; cross-halo terms are off. |
| Bulk velocity | **Withheld.** Each system is loaded at rest so it stays at its specified position. |
| Bulk position | Applied normally — the systems must be spatially separated. |
| Net momentum | Zeroed **per system** (this happens on every run, see below). |
| Halo monopole removal | Applied per system (halo centres are held static during warmup). |

**Video recording** is suppressed during warmup. Because the systems are isolated
and held at rest there, those frames are setup rather than simulation output, and
recording them would prepend a stretch of non-physical footage to every video. The
first frame written is the one at `t = 0`. For a 2.0-unit warmup at dt = 0.0005
that skips 4000 frames.

**At `t = 0`:** each system's stored bulk velocity is added to all of its
particles **and to its halo centre**, and full N-body coupling resumes. Adding a
uniform velocity leaves all internal relative motion untouched, so the relaxed
disc structure is preserved exactly — only the system as a whole starts moving.
From here each halo centre evolves as an inertial body under gravity rather than
being held in place (see `docs/dark-matter-halo.md`).

Because the bulk velocities are applied at `t = 0`, every orbital parameter in the
script still refers to `t = 0`. The warmup is prepended, not inserted.

#### How isolation is enforced under Octree gravity

Spatial separation alone does **not** isolate systems in a Barnes-Hut tree, and it
is worth being explicit about why:

1. The root node's bounding box spans every body, so it always encloses all
   systems regardless of their separation.
2. The multipole acceptance test (`d*d <= theta_sq * dsq`) gets *easier* to
   satisfy as distance grows — a distant system is **more** likely to be accepted
   as a single lump, not less. Greater separation makes the leak more efficient.
3. `BHNode` carries no system id, so the tree cannot filter by system.

Instead, the warmup builds a **separate tree per system**, containing only that
system's bodies, and evaluates only those bodies against it. Cross-system force is
then identically zero by construction, with no dependence on separation or on
`BH_Opening_Theta`.

This was verified numerically. Two systems loaded at rest and held through a
2-unit warmup drift by the same amount whether separated by 600, 3000, or 20000
code units (~3 units in every case), whereas genuine leakage would fall off as
1/r² (37 → 1.5 → 0.03 units over that range). A single system run alone wanders by
an identical 2.21 units, confirming the residual is each system's own
centre-of-mass wander from internal N-body noise, not cross-system coupling.

**Cost:** the tree is rebuilt once per system per step during warmup. With a
handful of systems this is comparable to the single-tree cost, since each tree
holds proportionally fewer bodies.

#### Choosing a value

`2.0` is a reasonable default — roughly 1.5 disc dynamical times for a galaxy with
a 210 km/s rotation curve, enough to damp the initial transients.

Be aware of a countervailing effect: an isolated disc **heats** through two-body
relaxation, so a long warmup raises the effective Toomre Q and leaves the disc
*less* responsive to the tidal perturbation. If you are chasing strong spiral-arm
contrast, prefer shorter warmups and treat long ones with suspicion.

**Default:** 0 (disabled)

---

## Automatic Corrections

Two corrections are applied unconditionally and have no script parameter. They are
documented here because they affect initial conditions.

### Per-System Net Momentum Removal

The procedural generators (`GalaxyDisc`, `SphericalUniverse`) zero their own
system's net momentum as part of building it, before applying any bulk velocity.

These generators sample particle azimuths randomly and never correct the sum,
leaving a residual net momentum of order `v_c / sqrt(N)`. Nothing damps it, so it
carries the system off the origin at constant velocity — for a 100k-particle disc
at 220 km/s that is ~1 code unit per time unit, enough to drift out of frame.

The correction lives in the generators themselves, so it applies exactly to the
systems that need it. A system assembled from explicit `Body` commands never
invokes it, which is correct: there the net momentum is physical (real
solar-system ephemerides, where the Sun's recoil balances the planets) and zeroing
it would corrupt the trajectories.

Ordering inside the generator matters and is deliberate:

1. Place particles with their internal motion only — no bulk velocity yet.
2. Zero the net momentum, so only sampling noise is removed.
3. Apply the bulk velocity (or, under `InitializationTime`, withhold it to `t = 0`).

Because step 2 precedes step 3, the intended bulk motion is never touched.

Applied per system rather than globally: a global correction cancels only the
*total*, leaving each individual system with its own residual drift, and under
warmup the isolated systems would wander away from their specified separation
before `t = 0`.

Subtracting a uniform velocity from a system changes no internal relative motion,
so it changes no physics — only which inertial frame that system starts in.

---

### `RemoveHaloMonopole` — Momentum Correction During Warmup

```
RemoveHaloMonopole  <0 | 1>
```

Controls a momentum-conservation correction for the rigid analytic halos **during warmup only**. When enabled (`1`, the default), each isolated system subtracts its own halo's net force (divided by its mass) uniformly from its particles' accelerations each step.

Background: a rigid analytic halo has no inertia, so on its own it cannot obey Newton's third law — a particle is pulled toward the halo centre but nothing pulls back. For a perfectly axisymmetric disc the per-particle forces cancel in the sum, but any asymmetry leaves a net force that would drift the whole system; the dominant offender is the **m=1 lopsided mode** (a two-armed m=2 bar is symmetric under 180° rotation and largely cancels). Subtracting a **uniform** vector from every particle leaves every difference `a_i - a_j` unchanged, so internal structure is preserved exactly and only the spurious bulk drift is cancelled — as if the rigid halo recoiled with its system.

**Post-warmup this correction is not used.** Once the systems couple at `t = 0`, each halo centre becomes an inertial body integrated under gravity (see `docs/dark-matter-halo.md`), and the net force the halo exerts on the particles is applied back onto the halo centre as its reaction. Momentum is then conserved directly by that recoil, so no monopole subtraction is wanted — applying one as well would double-count the reaction. During warmup the halo centres are instead held static (systems isolated and at rest), so the per-system subtraction remains the right tool there.

`FDE` is deliberately excluded from the correction, since it is a real external force whose net contribution is meant to be nonzero.

**Default:** 1 (enabled; warmup only).

**Default:** 1 (enabled)

---

### `RecordVideo` — Video Recording

```
RecordVideo  <0 | 1>
```

When enabled (`1`), records the simulation output to an MP4 video file using
FFmpeg, named after the script and written to an `output/` directory alongside it.
Each rendered frame is captured and encoded at 30 FPS.

Frames are recorded only while the simulation is running — pausing pauses
recording. If `InitializationTime` is set, the warmup frames (`t < 0`) are also
skipped, so the video begins at `t = 0`.

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

### `Camera_lookAt_System` — Follow a System

```
Camera_lookAt_System  <system_index>
```

Retargets the camera's look-at point onto the given system's **central body**
every frame, so a moving galaxy stays centred in view. `system_index` is 0-based
and must be a valid index into `N_SystemBodies`.

Only the look-at point moves. The camera's spherical offset from it — `phi`,
`theta` and radius — is preserved exactly, so:

- the viewing angle and zoom are unchanged as the target moves
- the runtime controls (W/A/S/D to orbit, J/L to zoom) keep working, now relative
  to the moving target
- `Camera_Orbit` composes with it: the camera orbits the followed body

The camera position is translated by the same delta as the look-at point, rather
than recomputed from stored angles, so following introduces no drift and does not
fight user input.

When this command is absent the feature is inactive and the fixed
`Camera_lookAt` point is used, exactly as before.

If both `Camera_lookAt` and `Camera_lookAt_System` are given, the system follow
wins — the fixed point is overwritten on the first frame.

**Errors.** The script is rejected with a message if the index is missing,
negative, or `>= N_Systems`. The upper bound is checked after the whole script is
read, so `Camera_lookAt_System` may appear before `N_SystemBodies`.

**Examples:**
```
Camera_lookAt_System  0     # follow the primary galaxy
Camera_lookAt_System  1     # follow the companion
```

Useful for an interaction where one galaxy has significant bulk motion: with a
fixed look-at the companion drifts out of frame, whereas following keeps it
centred for the whole run.

Note the target is the system's central body specifically, not its centre of mass.
For a galaxy these stay close: the central body sits near the halo
centre, though the two can differ by a kpc or two once strong tidal debris shifts
the system's barycentre away from the core.

**Default:** inactive (uses `Camera_lookAt`)

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
GalaxyDisc  <system>  <posX> <posY> <posZ>  <velX> <velY> <velZ>  <normalX> <normalY> <normalZ>  <central_body_mass> <disc_mass> <outer_radius> <inner_radius> <disc_scale_length> <toomre_Q>  <halo_circular_velocity> <halo_core_radius> <halo_truncation_radius>  [<sigma_z_ratio>]
```

Generates a flattened disc of particles with approximately circular orbits, representing a spiral galaxy. The disc plane is defined by the normal vector; particles orbit counter-clockwise when viewed from the direction the normal points.

| Parameter | Description |
|---|---|
| `system` | System index (0-based) |
| `posX`, `posY`, `posZ` | Position of the galaxy center |
| `velX`, `velY`, `velZ` | Bulk velocity of the entire galaxy |
| `normalX`, `normalY`, `normalZ` | Disc normal vector (defines orientation of the disc plane). Does not need to be unit length. Particles orbit counter-clockwise when viewed from the direction this vector points. Use `(0, 1, 0)` for a disc in the x-z plane |
| `central_body_mass` | Mass of the central body (the galaxy core / bulge), in code units. May be `0` for a pure disc; for a bulgeless model set it to a small token mass so it acts only as a centre marker |
| `disc_mass` | Total mass of the disc particles, in code units (must be `> 0`). Split evenly among the `N−1` disc particles at load time, so it is independent of the particle count |
| `outer_radius` | Maximum radius of the disc |
| `inner_radius` | Minimum radius of the disc (creates a central hole) |
| `disc_scale_length` | Exponential scale length of the disc, `h_r`. **Required**, and must satisfy `0 < h_r < outer_radius`. See below |
| `toomre_Q` | Target Toomre stability parameter. Controls the radial and tangential velocity dispersion via the Jeans equations. `Q = 1.0` is the stability threshold (disc will fragment); `Q = 1.2` gives a responsive disc with strong spiral structure; `Q = 1.5` gives a stable disc that responds only to external tidal perturbations; `Q = 2.0+` gives a hot, stable disc with weak spiral response. See `docs/toomre-q-velocity-dispersion.md` for the full derivation |
| `halo_circular_velocity` | Asymptotic circular velocity of the dark matter halo. Controls how flat the rotation curve is at large radii. Set to `0.0` for no halo |
| `halo_core_radius` | Core radius of the dark matter halo (cored isothermal sphere). The halo density flattens inside this radius. Irrelevant if `halo_circular_velocity` is 0 |
| `halo_truncation_radius` | Radius beyond which the halo's enclosed mass is frozen and the force falls off as `1/r^2`. `0` = untruncated. See below |
| `sigma_z_ratio` | **Optional** (20th field). Ratio of vertical to radial velocity dispersion, `sigma_z/sigma_r`, which sets the disc thickness via a self-gravitating sech² isothermal sheet and the vertical velocity spread. Defaults to `0.7` (Salo & Laurikainen 2000) when omitted, so existing 19-column scripts are unchanged. Solar-neighbourhood discs are nearer `0.5`; `0` gives a razor-thin, vertically cold disc |

#### Radial profile: scale length vs outer radius

Two radii describe the disc and they are independent:

| Quantity | Role |
|---|---|
| `disc_scale_length` (h_r) | How fast surface density falls: `Sigma(r) = Sigma_0 * exp(-r / h_r)`. Sets central concentration. |
| `outer_radius` (R) | Hard truncation. No particles are placed beyond it. |

`h_r` is the **e-folding distance** of the surface density: at `r = h_r` density is 37% of central, at `2 h_r` it is 13.5%, at `3 h_r` it is 5%. It sets not just appearance but the rotation curve shape (via enclosed mass) and the velocity dispersion (via the Toomre criterion).

A disc truncated at four scale lengths retains ~91% of its mass, which makes `R = 4 h_r` a common convention — but it is not universal, so both are required inputs. Salo & Laurikainen (2000) truncate M51a at `4 h_r` and M51b at `7.3 h_r`.

Particle radii are drawn from the exponential disc profile. The radial *number* density is the surface density times the area of a ring, `p(r) ~ r * exp(-r / h_r)`, which is a Gamma(shape 2, scale h_r) distribution. It is sampled exactly as `r = -h_r * ln(u1 * u2)` for two uniform deviates, with draws outside `[inner_radius, outer_radius]` rejected (~10% for a disc truncated at 4 h_r).

#### Halo truncation

`halo_truncation_radius` (Rh) sets where the halo's enclosed mass stops growing:

```
r <= Rh :  a_halo = Vc^2 * r / (r^2 + Rc^2)          cored isothermal
r >  Rh :  a_halo = M_halo(Rh) / r^2                 enclosed mass frozen
           M_halo(Rh) = Vc^2 * Rh^3 / (Rh^2 + Rc^2)
```

The two branches agree at `r = Rh`, so the force is continuous.

With `Rh = 0` the halo is untruncated and `M_halo(r) ~ Vc^2 * r` grows without bound at all radii. That overestimates long-range attraction, and in a two-galaxy encounter it distorts the effective mass ratio: each halo keeps accreting mass past its own disc edge, and the less concentrated halo gains proportionally more. For M51 the enclosed-mass ratio at the encounter separation comes out near 1:1.22 untruncated, versus the intended 1:1.82.

Setting `Rh` equal to the galaxy's own `outer_radius` makes the enclosed mass beyond the disc exactly the total mass within it, which is the prescription in Salo & Laurikainen (2000) section 2.2 (`Rh` = 400 arcsec for M51a, equal to its `Rd`).

The velocity dispersion at each radius is computed from the Toomre criterion: `sigma_r = Q * 3.36 * G * Sigma(r) / kappa(r)`, where `Sigma(r)` is the local surface density (exponential disc, normalized to the truncated mass) and `kappa(r)` is the epicyclic frequency derived from the full rotation curve (baryons + halo). The tangential dispersion follows from epicyclic theory: `sigma_phi = sigma_r * kappa / (2*Omega)`, and the mean streaming velocity is reduced below the circular speed by the asymmetric-drift correction. The vertical dispersion is `sigma_z = sigma_z_ratio * sigma_r` (default ratio 0.7), and particle heights follow a sech² isothermal sheet of scale `z0 = sigma_z^2 / (2*pi*G*Sigma)`. Together this produces a self-consistent equilibrium that suppresses particle-noise-driven instabilities while allowing the desired level of spiral response.

The dark matter halo applies a cored isothermal sphere potential: `a_halo = v_c^2 * r / (r^2 + r_c^2)` (optionally truncated beyond `halo_truncation_radius`, see above), centred on an **inertial halo centre** — a dynamical point integrated under gravity, carrying the halo's mass, rather than re-derived from the particles each step. Every particle feels its own system's halo plus the halo of every other system, so multiple galaxies interact through their halos as well as their particles; the reaction of that force is applied back to each halo centre (the "disc back-action"), so momentum is conserved and tidal debris does not drag the halo off the galaxy core. The halo is a rigid analytic background — its spherical shape never deforms and it is never tidally stripped, so it supplies no Chandrasekhar friction; see `docs/dark-matter-halo.md` for what this does and does not model.

**Example (Milky Way, disc in x-z plane).** `h_r` = 43.3 is the measured 2.6 kpc scale length, so `R / h_r` = 10.3:
```
GalaxyDisc  0   0.0 0.0 0.0   0.0 0.0 0.0   0.0 1.0 0.0   1500000.0 4500000.0 446.7 5.0 43.3 1.2  220.0 166.7 0.0
```

The same disc with an explicit, cooler vertical structure (`sigma_z/sigma_r` = 0.5, closer to the solar neighbourhood) appends the optional 20th field:
```
GalaxyDisc  0   0.0 0.0 0.0   0.0 0.0 0.0   0.0 1.0 0.0   1500000.0 4500000.0 446.7 5.0 43.3 1.2  220.0 166.7 0.0  0.5
```

**Example (M51b), a disc truncated at 7.3 scale lengths:**
```
GalaxyDisc  1   0.0 458.4 80.8   200.0 0.0 0.0   0.5373 0.8434 0.0000   780000.0 1170000.0 186.2 3.3 25.6 1.5  186.6 31.0 186.2
```

---

### `SphericalUniverse` — Procedural Spherical Distribution

```
SphericalUniverse  <system>  <posX> <posY> <posZ>  <velX> <velY> <velZ>  <total_mass> <radius> <H>  <halo_circular_velocity> <halo_core_radius> <halo_truncation_radius>
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
| `halo_truncation_radius` | Halo truncation radius; `0` for untruncated. Same meaning as for `GalaxyDisc` |

**Example (100k-body expanding sphere with Hubble flow for structure formation):**
```
SphericalUniverse  0   0.0 0.0 0.0   0.0 0.0 0.0   5.0e7 200.0 2.83  0.0 1.0 0.0
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

GalaxyDisc  0   0.0 0.0 0.0       0.0 0.0 0.0       0.0 1.0 0.0    1.0e7 5.0e6 250.0 25.0 62.5 1.2  200.0 50.0 0.0
GalaxyDisc  1   300.0 0.0 -500.0  -150.0 0.0 0.0    0.0 1.0 0.0    6.0e6 3.0e6 125.0 25.0 31.3 1.2  155.0 25.0 0.0
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
