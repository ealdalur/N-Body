# Toomre Q and Velocity Dispersion from the Jeans Equations

## 1. The Problem

A galaxy disc in N-body simulations is initialized with particles on circular orbits. But if the velocities are exactly circular, the disc is dynamically cold — any small perturbation (from particle noise, or a companion galaxy) causes the disc to fragment or develop spurious spiral structure that has nothing to do with the intended physics.

Real galaxy discs have **velocity dispersion** (random motions superimposed on the circular flow) that stabilizes them against gravitational collapse of local density perturbations. The Toomre Q parameter quantifies exactly how much dispersion is needed.

---

## 2. Toomre Stability Criterion

### 2.1 Derivation

Consider a patch of a razor-thin, self-gravitating disc with surface density Sigma, rotating with angular velocity Omega. A local density perturbation of wavelength lambda wants to collapse under its own gravity but is resisted by two effects:

1. **Rotation** (via the Coriolis force): perturbations larger than the epicyclic scale lambda_J = 4*pi^2*G*Sigma / kappa^2 are stabilized by differential rotation.

2. **Velocity dispersion** (pressure): perturbations smaller than the Jeans length lambda_J = sigma_r^2 / (G*Sigma) are stabilized by random motions.

If both scales overlap, ALL perturbations are stabilized. Toomre (1964) showed this requires:

```
Q = sigma_r * kappa / (3.36 * G * Sigma) >= 1
```

Where:
- `sigma_r` = radial velocity dispersion (km/s)
- `kappa` = epicyclic frequency (rad / time unit)
- `G` = gravitational constant
- `Sigma` = local surface density (mass / length^2)

### 2.2 Physical Interpretation

| Q value | Disc state |
|---------|-----------|
| Q < 1.0 | Unstable — disc fragments into bound clumps |
| Q = 1.0 | Marginally stable — any perturbation triggers collapse |
| Q = 1.2 | Cool disc — strong swing amplification, spontaneous multi-arm spirals |
| Q = 1.5 | Warm disc — stable to noise, responsive to external tidal perturbations |
| Q = 2.0 | Hot disc — only very strong perturbations produce visible structure |
| Q > 2.5 | Very hot — essentially featureless, no spiral response |

### 2.3 Solving for sigma_r

Rearranging the Toomre criterion:

```
sigma_r(r) = Q * 3.36 * G * Sigma(r) / kappa(r)
```

This gives the radial velocity dispersion at each radius that produces the target Q. The key insight: sigma_r is **not constant** across the disc — it varies with the local surface density and the local epicyclic frequency.

---

## 3. Epicyclic Frequency

### 3.1 Definition

The epicyclic frequency kappa describes how quickly a star oscillates radially about its guiding-center circular orbit. For a general rotation curve v_c(r):

```
kappa^2 = (2 * Omega / r) * d(r^2 * Omega) / dr
```

where Omega = v_c / r is the angular velocity.

### 3.2 Computation from v_c^2(r)

Since we know v_c^2 rather than Omega directly, it's more convenient to use:

```
kappa^2 = (1/r) * d(v_c^2)/dr + 2 * v_c^2 / r^2
```

This follows from expanding the definition with Omega^2 = v_c^2 / r^2.

### 3.3 For This Code's Rotation Curve

The circular velocity squared is:

```
v_c^2(r) = G * M_enc(r) / r + Vc^2 * r^2 / (r^2 + Rc^2)
```

where M_enc(r) is the enclosed baryonic mass (bulge + disc interior to r), and the second term is the halo contribution.

Taking the derivative:

```
d(v_c^2)/dr = -G * M_enc / r^2 + G * (dM_enc/dr) / r + Vc^2 * 2*r*Rc^2 / (r^2 + Rc^2)^2
```

For an exponential disc with scale length h_r:

```
dM_enc/dr = M_disc * (r / h_r^2) * exp(-r/h_r) / [1 - (1 + R/h_r)*exp(-R/h_r)]
```

This is the local surface density ring contribution at radius r.

### 3.4 Special Cases

- **Flat rotation curve** (v_c = const): kappa = sqrt(2) * Omega (Oort limit)
- **Solid body** (v_c proportional to r): kappa = 2 * Omega
- **Keplerian** (v_c proportional to 1/sqrt(r)): kappa = Omega

Most galaxy rotation curves are between flat and solid-body in the inner disc, transitioning to flat in the outer disc. So kappa varies between sqrt(2)*Omega and 2*Omega across the disc.

---

## 4. Tangential (Azimuthal) Velocity Dispersion

### 4.1 Epicyclic Theory

A star on a nearly-circular orbit traces an epicycle — an elliptical oscillation about the guiding center. In the rotating frame, the ratio of tangential to radial excursion is set by the Coriolis force:

```
sigma_phi / sigma_r = kappa / (2 * Omega)
```

This ratio varies across the disc:
- Flat rotation curve: sigma_phi/sigma_r = 1/sqrt(2) ~ 0.71
- Solid body: sigma_phi/sigma_r = 1.0
- Keplerian: sigma_phi/sigma_r = 0.5

### 4.2 Physical Meaning

The tangential dispersion is NOT independent of the radial dispersion. They are linked by the orbital mechanics of epicycles. You cannot set them independently without violating the Jeans equations (the disc would not be in equilibrium and would rearrange itself during the first few orbital periods).

This is why a flat Vtol (e.g., 10% perturbation to all velocity components equally) produces initial transients — it doesn't respect the epicyclic ratio and the disc has to relax into the correct sigma_phi/sigma_r before it settles.

---

## 5. Surface Density

### 5.1 Exponential Disc

The code uses an exponential surface density profile:

```
Sigma(r) = (M_disc / (2 * pi * h_r^2)) * exp(-r / h_r)
```

where h_r is the disc scale length (a required `GalaxyDisc` input, independent of
the truncation radius R) and M_disc = Mfrac * M is the total disc mass.

### 5.2 Why Sigma Matters

The Toomre criterion is local — Q depends on Sigma(r) at each radius. The exponential profile concentrates mass at small radii, so:
- Inner disc has high Sigma, high kappa → needs large sigma_r to maintain Q
- Outer disc has low Sigma, lower kappa → needs smaller sigma_r

This produces a physically realistic velocity dispersion profile that decreases outward, matching observations of real galaxies (where the inner disc is "hotter" than the outer disc).

---

## 6. The Jeans Equations (Context)

### 6.1 What They Are

The Jeans equations are the collisionless Boltzmann equation's first moments — they relate the density, mean velocity, and velocity dispersion tensor of a self-gravitating system. For a thin disc in steady state:

**Radial Jeans equation:**
```
d(nu * sigma_r^2)/dr + nu * (sigma_r^2 - sigma_phi^2) / r = -nu * dPhi/dr
```

where nu is the number density and Phi is the gravitational potential.

### 6.2 What This Code Does

Rather than solving the full Jeans equations iteratively (which requires knowing the distribution function), we use the **Toomre shortcut**:

1. Specify target Q (input parameter)
2. Compute sigma_r from Q, Sigma, kappa at each radius
3. Compute sigma_phi from the epicyclic ratio

This is exact for a thin disc and produces excellent equilibria. The full iterative Jeans approach (as in GalIC or AGAMA) additionally accounts for:
- Finite disc thickness (sech^2 vertical profile)
- Asymmetric drift correction (mean v_phi < v_c because pressure support partially replaces rotation)
- Higher-order moments

For 2D discs at moderate Q (1.0-2.0), the Toomre shortcut and the full Jeans solution agree to within a few percent.

### 6.3 Asymmetric Drift

When sigma_r is non-negligible, the mean tangential velocity is slightly less than v_c:

```
<v_phi> = v_c - sigma_r^2 / (2 * v_c) * [1 - (kappa/(2*Omega))^2 + d(ln Sigma)/d(ln r) + d(ln sigma_r^2)/d(ln r)]
```

This correction is typically 5-15 km/s for Q ~ 1.5 and is not currently implemented (particles are placed at v_c + Gaussian noise). For the purpose of tidal interaction simulations this is a second-order effect — the disc adjusts in the first ~0.5 dynamical times.

---

## 7. Implementation in This Codebase

### 7.1 Code Location

`Simulation.cpp`, function `LoadGalaxyDiscState()`.

### 7.2 Algorithm (per particle at radius r)

```cpp
// 1. Compute circular velocity (same as before)
double vc_sq = G * m_orbit / r + haloVc*haloVc * r*r / (r*r + haloRc_sq);
double vm = sqrt(vc_sq);
double Omega = vm / r;

// 2. Compute local surface density
double Sigma = (M_disc / (2*pi*h_r*h_r)) * exp(-r / h_r);

// 3. Compute epicyclic frequency from rotation curve derivative
double dMenc_dr = M_disc * (r / (h_r*h_r)) * exp(-r/h_r) / enc_denom;
double dvc_sq_dr = -G*m_orbit/(r*r) + G*dMenc_dr/r
                   + haloVc*haloVc * 2*r*haloRc_sq / ((r*r+haloRc_sq)*(r*r+haloRc_sq));
double kappa_sq = dvc_sq_dr/r + 2*vc_sq/(r*r);
double kappa = sqrt(kappa_sq);

// 4. Radial dispersion from Toomre Q
double sigma_r = Q * 3.36 * G * Sigma / kappa;

// 5. Tangential dispersion from epicyclic theory
double sigma_phi = sigma_r * kappa / (2*Omega);

// 6. Draw velocities from Gaussian distributions
double v_tan = -(vm + sigma_phi * normal(gen));   // tangential (mean = -vc, dispersion = sigma_phi)
double v_rad = sigma_r * normal(gen);              // radial (mean = 0, dispersion = sigma_r)
```

### 7.3 Sign Convention

- Tangential velocity is negative because `v_tan = -vm` gives counter-clockwise rotation when the tangential direction is defined as `cross(r_hat, n_hat)`.
- The tangential dispersion adds scatter around v_c, not around zero.
- The radial dispersion is centered on zero (no mean radial motion for a steady-state disc).

### 7.4 Edge Case: kappa_sq < 0

For very centrally concentrated mass distributions, the computed kappa_sq can go negative at extreme radii (numerical artifact from discrete enclosed mass). The code guards against this:

```cpp
if (kappa_sq < 0.0) kappa_sq = 4.0 * Omega * Omega;  // fall back to solid-body
```

This gives kappa = 2*Omega (solid body rotation), which is the maximum physically reasonable value and produces conservative (slightly high) velocity dispersion.

---

## 8. Tips for Setting Target Toomre Q

### 8.1 For Tidal Interaction Simulations (e.g., M51)

**Use Q = 1.5**

The galaxy that will be perturbed (M51a) needs to be:
- Stable enough that it doesn't develop spontaneous multi-arm structure from particle noise
- Responsive enough that the companion's tidal field produces visible two-arm spiral arms

Q = 1.5 hits this sweet spot. The m=2 tidal forcing is a coherent, strong perturbation — it works at any Q below ~2.0. Particle noise drives incoherent modes (m=2, 3, 4 all comparable) that require Q < ~1.2 to grow. So Q = 1.5 kills the noise while preserving the signal.

For the companion galaxy (M51b), Q = 1.5 is also fine — it keeps the companion stable during the encounter without affecting the tidal field it produces.

### 8.2 For Isolated Galaxy Evolution

**Use Q = 1.0-1.2**

If you want to study spontaneous spiral structure (swing amplification, bar formation, disc instabilities), use low Q:
- Q = 1.0: Marginally stable. Will develop strong multi-arm spirals and possibly a bar.
- Q = 1.2: Rich spiral structure from swing amplification. Multiple transient arms, possibly a weak bar.

### 8.3 For Galaxy Mergers / Collisions

**Use Q = 1.2-1.5**

Galaxy mergers involve violent tidal disruption. Q = 1.2 gives visually rich tidal tails and bridges. Q = 1.5 produces cleaner morphology (easier to see the tidal structure without background noise from swing amplification).

### 8.4 For Stability Tests

**Use Q = 2.0+**

If you want to verify that your disc remains axisymmetric (e.g., to test that the force solver isn't introducing artifacts), use Q = 2.0. At this value, the disc should remain featureless indefinitely in isolation.

### 8.5 Relationship to Particle Count

Unlike the old Vtol parameter, Toomre Q is **independent of particle count**. The physics is set by the continuous surface density and epicyclic frequency, not by the discrete particle properties. This means:
- Running the same simulation at 40k or 160k particles with the same Q gives the same macroscopic behavior
- You don't need to retune the stability parameter when changing resolution
- The softening length (r_soft) still needs to scale as 1/sqrt(N), but Q does not

### 8.6 Relationship to Softening

The effective Q is slightly boosted by gravitational softening (which suppresses the shortest-wavelength perturbations). For r_soft << h_r (which should always be the case), this correction is negligible and can be ignored.

### 8.7 Quick Reference

| Scenario | Recommended Q |
|----------|---------------|
| Tidal interaction (disc being perturbed) | 1.5 |
| Tidal interaction (perturber) | 1.5 |
| Isolated galaxy with spiral arms | 1.0-1.2 |
| Galaxy merger | 1.2-1.5 |
| Bar formation study | 1.0-1.1 |
| Stability / featureless disc | 2.0+ |
| "I see 3 arms before the encounter" | Increase Q to 1.5-1.8 |
| "The arms are too weak after encounter" | Decrease Q to 1.2-1.3 |
