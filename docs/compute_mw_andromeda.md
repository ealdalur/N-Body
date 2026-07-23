# compute_mw_andromeda.py

## Goal

Convert real astrophysical measurements of the Milky Way and Andromeda galaxies into the N-Body simulator's internal code unit system, producing ready-to-use parameters for the `Milky_Way_Andromeda.sim` script. This includes disc geometry, mass decomposition, dark matter halo parameters, relative kinematics, and Andromeda's disc orientation.

## Approach

1. Gather observational parameters for both galaxies from the literature.
2. Define the simulation unit system in terms of physical units.
3. Apply linear scaling to convert all dimensional quantities.
4. Compute Andromeda's disc normal vector from its observed inclination and position angle.
5. Estimate the collision timescale.
6. Output formatted GalaxyDisc script lines.

## Source Data

### Observational Parameters

| Parameter | Milky Way | Andromeda | Source |
|-----------|-----------|-----------|--------|
| Disc radius | 26.8 kpc | 33.5 kpc | Bland-Hawthorn & Gerhard 2016; Tamm et al. 2012 |
| Bulge mass | 1.5e10 Msun | 3.0e10 Msun | McMillan 2011; Tamm et al. 2012 |
| Disc mass | 4.5e10 Msun | 7.0e10 Msun | Bland-Hawthorn & Gerhard 2016; Tamm et al. 2012 |
| SMBH mass | 4.0e6 Msun (Sgr A*) | 1.4e8 Msun (M31*) | GRAVITY Collab. 2019; Bender et al. 2005 |
| Rotation curve | ~220 km/s (flat) | ~225 km/s (flat) | Eilers et al. 2019; Corbelli et al. 2010 |
| Halo core radius | 10 kpc | 12 kpc | Typical cored isothermal fits |
| Distance | 780 kpc | - | Riess et al. 2012 |
| Radial velocity | -110 km/s (approaching) | - | van der Marel et al. 2012 |
| Transverse velocity | 17 km/s | - | van der Marel et al. 2012 |
| Disc inclination | - | 77 deg | de Vaucouleurs 1958 |
| Position angle | - | 35 deg | de Vaucouleurs 1958 |

### Unit System

The simulation uses dimensionless "code units" defined by G=1:

| Dimension | 1 code unit = | Rationale |
|-----------|---------------|-----------|
| Distance | 60 pc = 0.060 kpc | Natural scale for galactic disc resolution with ~40k particles |
| Mass | 10^4 Msun | Sets disc particle masses to ~10^6 Msun (physically reasonable) |
| Velocity | 1 km/s | Direct read-off of rotation curves |
| Time | 58.7 Myr | Derived: du/vu = 60 pc / (1 km/s) = 58.7 Myr |

The simulation sets G=1 as a script parameter. With du=60 pc and vu=1 km/s fixed (distance from disc structure, velocity from rotation curves), tu=du/vu=58.7 Myr follows exactly. The mass unit is then determined by G=1: mu = du×vu²/G_SI ≈ 13,950 M☉. The commonly stated "10⁴ M☉" is a convenient approximation (accurate to ~18%).

## Equations and Calculations

### Unit Conversion

All conversions are simple linear scaling:

```
value_code = value_physical / unit_scale
```

For example:
```
mw_R = 26.8 kpc / 0.060 kpc = 446.7 code units
mw_M = 1.5e10 Msun / 1e4 Msun = 1,500,000 code units
```

### Mass Fraction (Mfrac)

The GalaxyDisc generator in the simulator creates N-1 disc particles whose total mass equals Mfrac times the central body mass M:

```
Mfrac = disc_mass / central_mass
MW: Mfrac = 4.5e10 / 1.5e10 = 3.0
And: Mfrac = 7.0e10 / 3.0e10 = 2.33
```

The central body mass represents the bulge (which subsumes the SMBH — the SMBH mass is negligible compared to the bulge).

### Andromeda Disc Normal Vector

The disc orientation is parameterized by inclination (tilt from line-of-sight) and position angle (rotation on the sky). The normal vector in the simulation frame is computed as:

```
n_x = cos(inclination)
n_y = sin(inclination) * cos(PA)
n_z = sin(inclination) * sin(PA)
```

With inclination = 77 deg, PA = 35 deg:
```
n_x = cos(77) = 0.2250
n_y = sin(77) * cos(35) = 0.7982
n_z = sin(77) * sin(35) = 0.5589
```

This produces a unit vector (|n| = 1 by construction) that tilts the disc away from the x-z plane, with the position angle rotating the tilt direction.

### Collision Timescale Estimate

A rough estimate assuming constant velocity:

```
t_collision ~ distance / v_radial
            = 780 kpc / (110 km/s * 1.022 kpc/Gyr per km/s)
            ~ 6.9 Gyr ~ 118 code time units
```

The factor 1.022 converts km/s to kpc/Gyr. This is a lower bound since gravitational acceleration increases the approach speed over time. The actual collision time from N-body simulations (e.g., van der Marel et al. 2012) is ~4.5 Gyr due to acceleration.

## Limitations and Assumptions

1. **Mass decomposition**: The bulge/disc split is observationally uncertain to ~30%. Different studies use different decomposition methods (photometric vs kinematic). The values used represent a consensus from multiple sources.

2. **Isothermal halo model**: The dark matter halo is modeled as a cored isothermal sphere with a_halo = v_c^2 * r / (r^2 + r_c^2). This is a simplification; real halos follow NFW profiles at large radii. The cored model avoids a central density cusp and matches flat rotation curves well within the disc region.

3. **Halo core radius**: The choice of r_c (10 kpc MW, 12 kpc And) is not well-constrained observationally. These values produce rotation curves that are approximately flat beyond a few kpc, which matches observations.

4. **Disc orientation simplification**: Andromeda's orientation is described by just two angles. The real disc has warps, rings, and non-axisymmetric structure that are not captured.

5. **Static initial conditions**: The galaxies are placed at their current observed positions with current observed velocities. No backward-integration to an earlier state is performed.

6. **Axisymmetric disc**: The GalaxyDisc generator produces a smooth, axisymmetric disc. Real galaxies have spiral arms, bars, and substructure.

7. **Particle resolution**: At 40,000 particles per galaxy, each disc particle represents ~10^6 Msun. Two-body relaxation is suppressed by the softening length but still exceeds the physical relaxation time by orders of magnitude.

## Code Implementation

### Structure

The script is organized into:
- **Unit definitions**: du, mu, vu establishing the code unit system.
- **Milky Way parameters**: Physical values in kpc, Msun, km/s.
- **Andromeda parameters**: Same set of physical values.
- **Relative geometry**: Distance, velocities, orientation angles.
- **Conversion and printing**: Each parameter is converted and printed alongside its physical equivalent.

### Unit Conversion Blocks

Each galaxy's parameters are converted in a straightforward block:

```python
mw_R = mw_radius_kpc / du        # kpc -> code distance
mw_M_central = mw_bulge_msun / mu  # Msun -> code mass
mw_haloVc = mw_vc_kms / vu        # km/s -> code velocity (trivially 1:1)
```

The velocity unit being 1 km/s means velocity values are numerically identical in physical and code units.

### Disc Normal Computation

The inclination/PA-to-normal conversion uses standard spherical-to-Cartesian mapping with the convention that inclination=0 means face-on to the observer (along x-axis in the sim frame). This produces the 3-component unit vector used by the GalaxyDisc generator to construct an orthonormal frame for particle placement.

### Collision Timescale

The rough timescale calculation uses the relation 1 km/s ~ 1.022 kpc/Gyr (from 1 pc = 3.086e13 km). This gives an order-of-magnitude estimate only; the actual time depends on gravitational acceleration during infall.

### Output

The final output section prints complete GalaxyDisc lines ready for copy-paste into a .sim file, including all 17 parameters in the correct order.
