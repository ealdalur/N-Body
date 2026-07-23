# compute_mw_andromeda_collision.py

## Goal

Determine the initial separation and approach velocity for a Milky Way - Andromeda collision scenario that begins interacting within a reasonable simulation runtime. The full-distance scenario (780 kpc apart, ~4.5 Gyr to collision) takes far too many timesteps. This script uses energy conservation to compute what the approach velocity would be at a much closer separation, as if the galaxies had already fallen inward from their current observed positions.

## Approach

1. Start from the observed boundary conditions: 780 kpc separation, 110 km/s radial approach velocity.
2. Estimate the total gravitating mass within a chosen close separation using the isothermal halo model.
3. Apply energy conservation (Keplerian two-body) to compute the velocity at the close separation.
4. Compare several candidate separations to choose one that gives a collision within a manageable number of timesteps while still allowing the galaxies to be initially non-overlapping.

## Source Data

### Initial Conditions (from observations)

| Parameter | Value | Source |
|-----------|-------|--------|
| Current separation | 780 kpc = 13,000 code units | Riess et al. 2012 |
| Current radial velocity | 110 km/s (approaching) | van der Marel et al. 2012 |
| MW disc radius | 446.7 code units (26.8 kpc) | Same as Milky_Way_Andromeda.sim |
| Andromeda disc radius | 558.3 code units (33.5 kpc) | Same as Milky_Way_Andromeda.sim |

### Galaxy Parameters (from `compute_mw_andromeda.py`)

| Parameter | MW | Andromeda |
|-----------|-----|-----------|
| Central mass (code) | 1,500,000 | 3,000,000 |
| Mfrac | 3.0 | 2.33 |
| Halo v_c | 220 km/s | 225 km/s |
| Halo r_c | 166.7 code | 200 code |

### Unit System

Same as `compute_mw_andromeda.py`:
- 1 distance unit = 60 pc = 0.060 kpc
- 1 mass unit = 10^4 Msun
- 1 velocity unit = 1 km/s
- 1 time unit = 58.7 Myr (= 60 pc / 1 km/s)

## Equations and Calculations

### Enclosed Mass Estimate

For a singular isothermal sphere (SIS), the enclosed mass within radius r is:

```
M_halo(r) = v_c^2 * r / G
```

With G=1 in code units, this simplifies to:

```
M_halo(r) = v_c^2 * r
```

At the chosen close separation r_close = 2500 (150 kpc), the total gravitating mass is estimated as:

```
M_total = M_halo_MW(r) + M_halo_And(r) + M_baryonic_MW + M_baryonic_And
```

where:
```
M_halo_MW = 220^2 * 2500 = 1.21e8 code units
M_halo_And = 225^2 * 2500 = 1.27e8 code units
M_bary_MW = 1,500,000 * (1 + 3.0) = 6,000,000 code units
M_bary_And = 3,000,000 * (1 + 2.33) = 9,990,000 code units
M_total ~ 2.56e8 code units
```

The baryonic contribution is only ~6% of the total — dark matter dominates the potential at these scales.

### Energy Conservation

Treating the two-galaxy system as a two-body problem, energy conservation gives:

```
(1/2)*v_0^2 - G*M/r_0 = (1/2)*v^2 - G*M/r
```

Solving for v at separation r (with G=1):

```
v^2 = v_0^2 + 2*M*(1/r - 1/r_0)
```

Substituting:
```
v^2 = 110^2 + 2 * 2.56e8 * (1/2500 - 1/13000)
v^2 = 12100 + 2 * 2.56e8 * (4e-4 - 7.7e-5)
v^2 = 12100 + ~165,000
v ~ 420 km/s
```

(The actual script uses slightly different mass estimates depending on the separation comparison table.)

### Time to First Contact

A naive estimate of how long until the disc edges overlap:

```
gap = separation - R_MW - R_And
t_contact ~ gap / v_approach
```

For r_close = 3000 (the value chosen for the script):
```
gap = 3000 - 446.7 - 558.3 = 1995 code units
t_contact ~ 1995 / 300 ~ 6.7 code time units ~ 393 Myr
```

In practice, gravitational acceleration shortens this further.

### Separation Comparison Table

The script evaluates multiple candidate separations to help choose the best trade-off:

| Separation | Gap to overlap | Approach velocity | Time to contact |
|------------|---------------|-------------------|-----------------|
| 5000 | 3995 | lower | longer |
| 4000 | 2995 | moderate | moderate |
| 3500 | 2495 | moderate | moderate |
| 3000 | 1995 | ~300 km/s | ~393 Myr |
| 2500 | 1495 | higher | shortest |

The final script uses r = 3000 (180 kpc) as a balance: the discs are well-separated initially (giving time for the halos to interact before disc overlap) but close enough that first passage occurs within ~10 simulation time units.

### Transverse Velocity

A transverse component of 50 km/s is added to produce a slightly off-center (grazing) encounter rather than a head-on collision. This is both more physically realistic (the observed proper motion implies some tangential velocity) and produces more visually interesting tidal tails and asymmetric merger dynamics.

## Limitations and Assumptions

1. **Point-mass energy conservation**: The two-body energy conservation formula treats both galaxies as point masses. In reality, the extended mass distribution means the effective gravitating mass changes as the galaxies approach each other. This makes the computed velocity an estimate (probably an overestimate since less mass is actually enclosed at the midpoint between the galaxies).

2. **Singular isothermal sphere**: The enclosed mass formula M = v_c^2 * r assumes an SIS profile extending to infinity. The simulation uses a cored isothermal profile (M ~ v_c^2 * r only for r >> r_c), which reduces the enclosed mass at small radii. At r=2500-3000 >> r_c=167-200, this is a good approximation.

3. **No dynamical friction estimate**: As the galaxies approach, dynamical friction from overlapping halos decelerates them. The script does not model this — the simulation's cross-halo gravitational terms handle it self-consistently during runtime.

4. **Separation comparison uses half the halo mass**: The comparison table uses `(v_c_MW^2 + v_c_And^2) * r * 0.5` as an effective mass estimate. The factor of 0.5 is a rough correction for the fact that at separation r, each galaxy's halo extends only to r from its own center (not the other's center). This is a crude approximation.

5. **Instantaneous placement**: The galaxies are placed at their positions with the computed velocities instantaneously. No integration history is preserved — the disc structure is that of a relaxed, equilibrium galaxy, not one that has been tidally pre-deformed during infall. This is reasonable since tidal effects are weak until the galaxies are within a few disc radii of each other.

6. **Fixed transverse velocity**: The 50 km/s transverse component is a free parameter chosen for visual quality, not derived from observations. The real MW-Andromeda encounter likely has a smaller impact parameter.

## Code Implementation

### Structure

- **Parameter block**: Repeats the galaxy parameters from the main MW-Andromeda script (code units).
- **Mass estimation**: Computes enclosed halo mass at the target separation.
- **Energy conservation**: Solves for approach velocity.
- **Comparison table**: Loops over candidate separations for decision-making.
- **Output**: Prints the chosen parameters in script-ready format.

### Mass Estimation

```python
mw_halo_enclosed = mw_haloVc**2 * r_close
and_halo_enclosed = and_haloVc**2 * r_close
mw_baryonic = mw_M * (1 + mw_Mfrac)
and_baryonic = and_M * (1 + and_Mfrac)
M_total = mw_halo_enclosed + and_halo_enclosed + mw_baryonic + and_baryonic
```

The halo mass formula is the SIS enclosed mass with G=1. The baryonic mass includes both the central body (bulge) and disc particles (bulge * (1 + Mfrac) = total baryonic).

### Velocity Calculation

```python
v_sq = v0**2 + 2 * 1.0 * M_total * (1.0/r_close - 1.0/r0)
v_approach = math.sqrt(v_sq)
```

Direct application of energy conservation. The factor `1.0` is G (=1 in code units). Since 1/r_close > 1/r0 (closer = deeper potential), the term is positive, increasing the velocity from its initial value.

### Separation Comparison Loop

The loop tries multiple separations with a modified mass estimate:

```python
M_eff = (mw_haloVc**2 + and_haloVc**2) * r * 0.5 + mw_baryonic + and_baryonic
```

The 0.5 factor accounts for the geometric reality that at separation r, each galaxy's halo is not fully enclosed within r of the system barycenter. This gives a more conservative (lower) velocity estimate than using the full M_total, which is appropriate for the comparison since the goal is to avoid overestimating the approach speed.

### Output

The script prints both diagnostic information (mass breakdowns, energy terms) and final values suitable for constructing the .sim script. The camera is placed at the midpoint between galaxies with a y-offset for an elevated viewing angle.
