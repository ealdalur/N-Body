# Keplerian Orbital Mechanics: A First-Principles Derivation

## Introduction

Any two bodies interacting solely through Newtonian gravity trace conic sections (circles, ellipses, parabolas, or hyperbolas) about their common center of mass. For bound orbits (negative total energy), the trajectory is an ellipse. This document derives, from first principles, how six numbers — the Keplerian orbital elements — fully specify the shape, orientation, and phase of an elliptical orbit, and how to convert those six numbers into a Cartesian position and velocity vector at any instant.

---

## 1. The Two-Body Problem

### 1.1 Equations of Motion

Consider two point masses m_1 and m_2 separated by vector **r** = **r**_2 - **r**_1. Newton's law of gravitation gives:

```
m_1 * d²r_1/dt² = +G*m_1*m_2 / |r|² * r_hat
m_2 * d²r_2/dt² = -G*m_1*m_2 / |r|² * r_hat
```

Subtracting (dividing each by its respective mass) gives the relative acceleration:

```
d²r/dt² = -G*(m_1 + m_2) / |r|² * r_hat = -mu / r² * r_hat
```

where **mu = G*(m_1 + m_2)** is the gravitational parameter. This is a central force problem — the acceleration is always directed toward the origin (the other body).

### 1.2 Conservation Laws

Two conservation laws constrain the motion:

**Angular momentum:**
```
L = r × v = constant vector
```

Since **L** is constant, the motion is confined to a plane (the orbital plane, perpendicular to **L**). This reduces the 3D problem to 2D.

**Energy:**
```
E = (1/2)|v|² - mu/r = constant
```

For bound orbits, E < 0.

### 1.3 The Orbit Equation

Working in polar coordinates (r, f) within the orbital plane, where f is the true anomaly (angle from periapsis), the solution to the equations of motion is:

```
r = a(1 - e²) / (1 + e*cos(f))
```

This is the equation of a conic section with:
- **a** = semi-major axis (size of the ellipse)
- **e** = eccentricity (shape of the ellipse)
- **f** = true anomaly (angular position measured from periapsis)

For 0 <= e < 1, this is an ellipse. The closest approach (periapsis, f=0) has r = a(1-e), and the farthest point (apoapsis, f=pi) has r = a(1+e).

---

## 2. The Six Keplerian Elements

An elliptical orbit in 3D space has six degrees of freedom — matching the six components of an initial position and velocity vector. The Keplerian elements decompose these six numbers into physically meaningful quantities:

### 2.1 Shape Elements (define the ellipse geometry)

**Semi-major axis (a):**
The half-length of the ellipse's longest diameter. Determines the orbital energy and period:
```
E = -mu / (2a)
P = 2*pi * sqrt(a³ / mu)
```
Units: length (AU for solar system work).

**Eccentricity (e):**
Measures departure from circularity. For ellipses, 0 <= e < 1.
- e = 0: perfect circle (r = a everywhere)
- e = 0.5: moderate ellipse (apoapsis is 3x periapsis)
- e -> 1: increasingly elongated

The distance from center to focus is a*e. Periapsis distance is a*(1-e), apoapsis is a*(1+e).

### 2.2 Orientation Elements (define the orbital plane's attitude in space)

These three angles orient the orbital plane relative to a reference frame (the ecliptic for solar system work):

**Inclination (i):**
The angle between the orbital plane and the reference plane (ecliptic). Ranges from 0 to 180 deg.
- i = 0: orbit lies in the ecliptic, same direction as Earth
- i = 90: orbit is perpendicular to the ecliptic
- i > 90: retrograde orbit (opposite to Earth's motion)

**Longitude of the ascending node (Omega):**
The angle, measured in the reference plane from the reference direction (vernal equinox), to where the orbit crosses the reference plane going "up" (from below to above the ecliptic). Ranges from 0 to 360 deg.

The ascending node is the point where the body crosses from negative z (below ecliptic) to positive z (above ecliptic). Omega specifies which direction this crossing point lies.

**Argument of periapsis (omega):**
The angle, measured in the orbital plane from the ascending node to the periapsis point, in the direction of orbital motion. Ranges from 0 to 360 deg.

This locates the closest-approach point within the tilted orbital plane.

### 2.3 Phase Element (locates the body along the orbit at a specific time)

**Mean anomaly (M):**
A fictitious angle that advances uniformly with time:
```
M(t) = M_0 + n*(t - t_0)
```
where n = 2*pi/P is the mean motion (constant angular rate).

M is not the actual angular position of the body (that is the true anomaly f). Instead, M parameterizes time in angular units: M = 0 at periapsis passage, M = pi at apoapsis, M = 2*pi after one full period.

The mean anomaly advances at a constant rate regardless of where the body is in its orbit, making it the natural time-like variable for orbital computations.

---

## 3. The Solar System Reference Frame

### 3.1 The Ecliptic Coordinate System

NASA, JPL, and the IAU use the **International Celestial Reference Frame (ICRF)** for high-precision work, but for orbital element catalogs, the standard frame is the **ecliptic and mean equinox of J2000.0**:

**Reference epoch:** J2000.0 = January 1, 2000, 12:00 TT (Julian Date 2451545.0)

**Reference plane:** The ecliptic — the plane of Earth's orbit around the Sun. This is the natural reference for solar system dynamics because most planets orbit close to this plane.

**Reference direction:** The vernal equinox (First Point of Aries) — the direction from the Sun to Earth at the March equinox. In the ecliptic frame, this defines the x-axis.

**Axes:**
```
x: toward vernal equinox (Aries direction)
y: in the ecliptic plane, 90 deg ahead of x in Earth's orbital direction
z: toward the north ecliptic pole (perpendicular to ecliptic, right-hand rule)
```

### 3.2 Why J2000.0?

The ecliptic plane and equinox direction slowly drift due to precession (Earth's spin axis traces a 26,000-year cone) and nutation. Fixing the frame at J2000.0 provides a non-rotating reference that all measurements can be reduced to, regardless of when the observation was made.

### 3.3 Usage in Space Missions

JPL's HORIZONS system reports state vectors in the ecliptic J2000 frame. NASA mission planning (trajectory design, navigation) uses this frame for interplanetary work. The SPICE toolkit provides transformations between J2000 ecliptic, J2000 equatorial (ICRF), and body-fixed frames.

Orbital elements published for planets (such as those from Standish 1992, used by this simulator) are referenced to the J2000 ecliptic, with angles measured from the vernal equinox.

---

## 4. Kepler's Equation

### 4.1 The Problem

Given the mean anomaly M (which we can compute at any time from M = M_0 + n*t), we need to find the body's actual position. The connection requires an intermediate variable: the **eccentric anomaly E**.

### 4.2 Geometric Derivation

Consider the orbit ellipse with semi-major axis a and semi-minor axis b = a*sqrt(1-e²). Circumscribe a circle of radius a around the same center. For any point P on the ellipse at angle f (true anomaly) from periapsis:

1. Draw a perpendicular from P to the major axis (the x-axis of the orbital plane).
2. Extend this perpendicular up to where it intersects the circumscribed circle at point Q.
3. The angle from the center of the ellipse to Q, measured from the periapsis direction, is the eccentric anomaly E.

The geometric relationship between the point on the ellipse and the corresponding point on the circle gives:
```
x_ellipse = a*cos(E) - a*e = a*(cos(E) - e)
y_ellipse = b*sin(E) = a*sqrt(1-e²)*sin(E)
```

(This will be derived in detail in Section 5.)

### 4.3 Area Relationship

Kepler's second law states that equal areas are swept in equal times. The rate of area swept is:
```
dA/dt = L / (2*m) = constant
```

For the full orbit:
```
A_ellipse = pi*a*b = (L / (2*m)) * P
```

Now consider the area swept from periapsis to the current position. In the circumscribed circle, the corresponding area (from center to Q through angle E) is:

```
A_circle_sector = (1/2)*a²*E                    [sector from 0 to E]
A_triangle = (1/2)*a²*e*sin(E)                  [triangle from focus to vertical]
A_circle(from focus) = (1/2)*a²*(E - e*sin(E))  [sector minus triangle]
```

The area on the ellipse relates to the area on the circle by the ratio b/a (since the ellipse is a vertically compressed circle):
```
A_ellipse(from focus) = (b/a) * (1/2)*a²*(E - e*sin(E))
                      = (1/2)*a*b*(E - e*sin(E))
```

The fraction of the period elapsed is:
```
t/P = A_swept / A_total = (1/2)*a*b*(E - e*sin(E)) / (pi*a*b) = (E - e*sin(E)) / (2*pi)
```

Since M = 2*pi*t/P by definition:
```
M = E - e*sin(E)
```

This is **Kepler's equation**. It relates the time-like variable M to the geometric variable E.

### 4.4 Interpretation

- At periapsis: E = 0, M = 0 (sin(0) = 0, consistent)
- At apoapsis: E = pi, M = pi (sin(pi) = 0, consistent)
- In between: E leads M (the body moves faster near periapsis than mean motion would suggest)

For e = 0 (circular orbit): M = E = f (all three anomalies are identical and advance uniformly).

### 4.5 Solving Kepler's Equation

Kepler's equation is transcendental — E appears both as itself and inside sin(E). There is no closed-form solution. The standard numerical approach is Newton-Raphson iteration:

Define: g(E) = E - e*sin(E) - M = 0

Newton's method:
```
E_{n+1} = E_n - g(E_n) / g'(E_n)
         = E_n - (E_n - e*sin(E_n) - M) / (1 - e*cos(E_n))
```

Starting from E_0 = M (which is exact for e=0 and a good approximation for small e):
- For e < 0.3 (all planets): converges in 3-5 iterations to machine precision
- For e near 1 (comets): may need ~10 iterations or a better starting guess

The denominator (1 - e*cos(E)) is always positive for e < 1, so the iteration is unconditionally stable for elliptical orbits.

---

## 5. Position in the Orbital Plane

### 5.1 Derivation from the Eccentric Anomaly

Once E is known, the position in the orbital plane (with the x-axis pointing from the focus toward periapsis) follows directly from the circumscribed-circle construction.

**x-coordinate:** The x-position on the circumscribed circle is a*cos(E), measured from the ellipse center. The focus is offset from the center by a*e (toward periapsis). So the distance from the focus, measured along x, is:

```
x = a*cos(E) - a*e = a*(cos(E) - e)
```

**y-coordinate:** The y-position on the circle is a*sin(E). The ellipse is a circle compressed by factor b/a = sqrt(1-e²) in the y-direction:

```
y = a*sqrt(1-e²) * sin(E)
```

### 5.2 Verification: Distance Formula

The radial distance from the focus should match the orbit equation r = a(1-e²)/(1+e*cos(f)). Compute:

```
r² = x² + y²
   = a²*(cos(E) - e)² + a²*(1-e²)*sin²(E)
   = a²*(cos²(E) - 2e*cos(E) + e² + sin²(E) - e²*sin²(E))
   = a²*(1 - 2e*cos(E) + e² - e²*sin²(E))
   = a²*(1 - 2e*cos(E) + e²*cos²(E))
   = a²*(1 - e*cos(E))²
```

Therefore:
```
r = a*(1 - e*cos(E))
```

This is the standard radial distance formula in terms of E, confirming our position expressions are correct.

### 5.3 Relationship Between E and f

The true anomaly f (actual angular position) relates to E by:
```
tan(f/2) = sqrt((1+e)/(1-e)) * tan(E/2)
```

This can be derived by expressing cos(f) in terms of x/r:
```
cos(f) = x/r = (cos(E) - e) / (1 - e*cos(E))
```

For the simulation code, we do not need f explicitly — we work directly with the (x, y) Cartesian coordinates.

---

## 6. Velocity in the Orbital Plane

### 6.1 Derivation from Time Derivatives

The velocity requires differentiating position with respect to time. We use the chain rule through E:

```
dx/dt = (dx/dE) * (dE/dt)
dy/dt = (dy/dE) * (dE/dt)
```

**Spatial derivatives:**
```
dx/dE = -a*sin(E)
dy/dE = a*sqrt(1-e²)*cos(E)
```

**Time derivative of E:** Differentiate Kepler's equation M = E - e*sin(E) with respect to time:
```
dM/dt = dE/dt - e*cos(E)*dE/dt = dE/dt * (1 - e*cos(E))
```

Since dM/dt = n (the mean motion, a constant):
```
dE/dt = n / (1 - e*cos(E))
```

### 6.2 Velocity Components

Combining:
```
vx = dx/dt = -a*sin(E) * n / (1 - e*cos(E))
           = -a*n*sin(E) / (1 - e*cos(E))

vy = dy/dt = a*sqrt(1-e²)*cos(E) * n / (1 - e*cos(E))
           = a*n*sqrt(1-e²)*cos(E) / (1 - e*cos(E))
```

where n = sqrt(mu/a³) is the mean motion.

### 6.3 Verification: Speed at Periapsis and Apoapsis

At periapsis (E=0):
```
vx = 0 (tangent to orbit, no radial component)
vy = a*n*sqrt(1-e²) / (1-e)
```

The vis-viva equation gives v² = mu*(2/r - 1/a). At periapsis, r = a*(1-e):
```
v² = mu*(2/(a*(1-e)) - 1/a) = mu/(a) * (2/(1-e) - 1) = mu/(a) * (1+e)/(1-e)
```

Our formula gives:
```
vy² = a²*n²*(1-e²)/(1-e)² = (mu/a)*(1-e²)/(1-e)² = (mu/a)*(1+e)/(1-e)
```

These match, confirming the velocity derivation.

At apoapsis (E=pi):
```
vx = 0
vy = -a*n*sqrt(1-e²) / (1+e)
```

The speed is lower at apoapsis (Kepler's second law: slower when farther from the Sun).

---

## 7. Rotation to the Ecliptic Frame

### 7.1 The Problem

The position (x, y) and velocity (vx, vy) computed above are in the orbital plane, with:
- x-axis: toward periapsis (from the focus)
- y-axis: perpendicular, in the direction of motion at periapsis
- z-axis: along the angular momentum vector (perpendicular to the orbital plane)

We need to rotate this into the ecliptic frame, where:
- x-axis: toward vernal equinox
- y-axis: 90 deg ahead in the ecliptic
- z-axis: toward north ecliptic pole

This requires three rotations, defined by the three orientation angles (omega, i, Omega).

### 7.2 The Three Rotations

The rotation from orbital-plane coordinates to ecliptic coordinates is decomposed into three steps, applied in order:

**Step 1: Rotate by omega (argument of periapsis) about the z-axis of the orbital plane.**

This rotates from the periapsis direction to the ascending node direction. After this rotation, the x-axis points toward the ascending node.

```
R_z(omega) = | cos(omega)  -sin(omega)  0 |
             | sin(omega)   cos(omega)  0 |
             |     0            0       1 |
```

**Step 2: Rotate by i (inclination) about the new x-axis (the node line).**

This tilts the orbital plane relative to the ecliptic. After this rotation, the z-axis is aligned with the ecliptic pole.

```
R_x(i) = | 1     0       0    |
          | 0   cos(i)  -sin(i)|
          | 0   sin(i)   cos(i)|
```

**Step 3: Rotate by Omega (longitude of ascending node) about the new z-axis (ecliptic pole).**

This rotates the node line from the x-axis to its actual ecliptic longitude.

```
R_z(Omega) = | cos(Omega)  -sin(Omega)  0 |
             | sin(Omega)   cos(Omega)  0 |
             |     0            0        1 |
```

### 7.3 Combined Rotation Matrix

The full transformation is:

```
r_ecliptic = R_z(Omega) * R_x(i) * R_z(omega) * r_orbital
```

Since r_orbital = (x', y', 0) (the motion is in-plane, z'=0), we only need the first two columns of the combined matrix. Multiplying out:

```
R = R_z(Omega) * R_x(i) * R_z(omega)
```

Let: cos_o = cos(omega), sin_o = sin(omega), cos_O = cos(Omega), sin_O = sin(Omega), cos_i = cos(i), sin_i = sin(i).

The product (keeping only terms that multiply x' and y', since z'=0):

**Column 1 (multiplies x'):**
```
R_11 = cos_O*cos_o - sin_O*sin_o*cos_i
R_21 = sin_O*cos_o + cos_O*sin_o*cos_i
R_31 = sin_o*sin_i
```

**Column 2 (multiplies y'):**
```
R_12 = -cos_O*sin_o - sin_O*cos_o*cos_i
R_22 = -sin_O*sin_o + cos_O*cos_o*cos_i
R_32 = cos_o*sin_i
```

### 7.4 Final Ecliptic Coordinates

```
x_ecl = R_11*x' + R_12*y'
y_ecl = R_21*x' + R_22*y'
z_ecl = R_31*x' + R_32*y'
```

The same matrix applies to the velocity:
```
vx_ecl = R_11*vx' + R_12*vy'
vy_ecl = R_21*vx' + R_22*vy'
vz_ecl = R_31*vx' + R_32*vy'
```

This works because the rotation matrix is constant (the orbital elements don't change with time in a two-body problem), so d(R*r)/dt = R*dr/dt.

### 7.5 Derivation of Matrix Elements (Example)

To show how the matrix elements arise, consider R_11:

```
R_z(Omega) * R_x(i) * R_z(omega) applied to (1, 0, 0)^T:
```

First, R_z(omega) * (1,0,0)^T = (cos_o, sin_o, 0)^T

Then, R_x(i) * (cos_o, sin_o, 0)^T = (cos_o, sin_o*cos_i, sin_o*sin_i)^T

Finally, R_z(Omega) * (cos_o, sin_o*cos_i, sin_o*sin_i)^T:
```
x = cos_O*cos_o - sin_O*sin_o*cos_i
y = sin_O*cos_o + cos_O*sin_o*cos_i
z = sin_o*sin_i
```

This gives R_11, R_21, R_31. The second column follows identically with initial vector (0, 1, 0)^T.

### 7.6 Special Cases

**i = 0 (equatorial orbit):** cos_i = 1, sin_i = 0. The z-component vanishes and:
```
x_ecl = cos(Omega + omega)*x' - sin(Omega + omega)*y'
y_ecl = sin(Omega + omega)*x' + cos(Omega + omega)*y'
z_ecl = 0
```
Only the sum (Omega + omega) matters — individual values are degenerate (there is no ascending node for an orbit in the reference plane).

**omega = 0 (periapsis at ascending node):**
```
R_11 = cos_O
R_21 = sin_O
R_31 = 0
```
The x-axis of the orbital plane aligns with the node line.

---

## 8. Frame Conversion: Z-Pole Ecliptic to Y-Pole Ecliptic

### 8.1 Motivation

The standard ecliptic frame (used by NASA, JPL, all orbital element catalogs) has:
```
x: toward vernal equinox
y: 90 deg ahead in the ecliptic
z: north ecliptic pole
```

The N-Body simulator uses a different convention:
```
x: toward vernal equinox (same)
y: north ecliptic pole (UP in the display)
z: 90 deg BEHIND in the ecliptic (so orbits go CCW in x-z viewed from +y)
```

The simulator's choice makes the ecliptic plane coincide with the x-z horizontal plane, with y as the vertical/display axis. This is a common convention in 3D graphics where y is "up."

### 8.2 Requirements for the Transformation

We need a rotation (proper, det = +1) that maps:
```
Standard basis:              Simulator basis:
x_std = (1,0,0)     ->      x_sim = (1,0,0)     [x stays the same]
y_std = (0,1,0)     ->      z_sim = (0,0,-1)    [ecliptic y maps to -z in sim]
z_std = (0,0,1)     ->      y_sim = (0,1,0)     [ecliptic pole maps to +y in sim]
```

### 8.3 Derivation as a Rotation About the X-Axis

Since x is unchanged, this must be a rotation about the x-axis. The general form of R_x(theta) is:

```
R_x(theta) = | 1    0          0      |
             | 0    cos(theta) -sin(theta) |
             | 0    sin(theta)  cos(theta) |
```

We need:
```
y_std = (0, 1, 0) -> (0, cos(theta), sin(theta)) = (0, 0, -1)
```

This requires cos(theta) = 0 and sin(theta) = -1, giving **theta = -90 deg = -pi/2**.

Verification:
```
R_x(-pi/2) = | 1   0    0 |
             | 0   0    1 |
             | 0  -1    0 |
```

Apply to each standard basis vector:
```
(1, 0, 0) -> (1, 0, 0)   [x unchanged]       ✓
(0, 1, 0) -> (0, 0, -1)  [y -> -z]           ✓
(0, 0, 1) -> (0, 1, 0)   [z -> +y]           ✓
```

### 8.4 Verifying It Is a Proper Rotation

A proper rotation has det(R) = +1 (preserves orientation/handedness):

```
det(R_x(-pi/2)) = 1 * (0*0 - 1*(-1)) = 1 * 1 = +1   ✓
```

The transformation preserves:
- Handedness (right-hand rule)
- Distances (|Rv| = |v| for all v)
- Angles between vectors (v·w = (Rv)·(Rw))
- Cross products: R(a × b) = (Ra) × (Rb)

This means physical quantities like angular momentum direction and orbital sense (prograde vs retrograde) are preserved. An orbit that is counter-clockwise in the standard frame remains counter-clockwise in the simulator frame.

### 8.5 The Transformation in Component Form

For any vector in the standard ecliptic frame (x_ecl, y_ecl, z_ecl), the simulator frame coordinates are:

```
x_sim = x_ecl
y_sim = z_ecl
z_sim = -y_ecl
```

And the inverse (simulator to standard):
```
x_ecl = x_sim
y_ecl = -z_sim
z_ecl = y_sim
```

### 8.6 Effect on Orbital Motion

In the standard frame, a planet in the ecliptic plane (z=0) with prograde motion has:
- Position: (x, y, 0) tracing an ellipse in the x-y plane
- Velocity: (vx, vy, 0) with angular momentum pointing toward +z

After transformation:
- Position: (x, 0, -y) tracing an ellipse in the x-z plane
- Velocity: (vx, 0, -vy) with angular momentum pointing toward +y

Viewed from above (+y direction in the simulator), the orbit goes counter-clockwise in the x-z plane — the same sense as in the standard frame viewed from +z.

### 8.7 Why Not a Coordinate Relabeling?

One might be tempted to simply relabel axes (call z "y" and y "-z"). While this produces the same numbers, it is an **improper** transformation (det = -1, a reflection). Physically, this would flip the handedness of the coordinate system, reversing the sense of cross products and angular momentum.

Using R_x(-90 deg) instead ensures that:
- The cross product formula remains valid
- Angular velocities point in the correct direction
- The right-hand rule still applies
- Torques and angular momentum are consistent

This matters for the simulator because forces (like gravity) are computed using cross products (for angular momentum checks) and the leapfrog integrator assumes a right-handed coordinate system.

---

## 9. Complete Pipeline Summary

Given the six Keplerian elements (a, e, i, Omega, omega, M) and the gravitational parameter mu at a reference epoch:

1. **Solve Kepler's equation** M = E - e*sin(E) for E using Newton-Raphson.

2. **Position in orbital plane:**
   ```
   x' = a*(cos(E) - e)
   y' = a*sqrt(1-e²)*sin(E)
   ```

3. **Velocity in orbital plane:**
   ```
   n = sqrt(mu/a³)
   vx' = -a*n*sin(E) / (1 - e*cos(E))
   vy' = a*n*sqrt(1-e²)*cos(E) / (1 - e*cos(E))
   ```

4. **Rotate to ecliptic frame** using R_z(Omega)*R_x(i)*R_z(omega) applied to both (x', y', 0) and (vx', vy', 0).

5. **Rotate to simulator frame** using R_x(-90 deg):
   ```
   (x_sim, y_sim, z_sim) = (x_ecl, z_ecl, -y_ecl)
   ```

The result is a 6-component state vector (x, y, z, vx, vy, vz) in the simulator's coordinate system, fully specifying the body's position and velocity at the reference epoch.

---

## 10. Propagation to Other Times

To find the state at a later time t:

1. Advance the mean anomaly: M(t) = M_0 + n*(t - t_0)
2. Re-solve Kepler's equation for E(t)
3. Recompute position and velocity

All other elements (a, e, i, Omega, omega) remain constant in a pure two-body problem. In reality, planetary perturbations cause these elements to drift slowly (secular perturbations) — this is why the N-Body simulator propagates the full equations of motion rather than relying on Keplerian propagation.

---

## References

- Kepler, J. (1609). *Astronomia Nova*. (Introduction of the ellipse and area law.)
- Newton, I. (1687). *Principia Mathematica*. (Derivation of Kepler's laws from gravity.)
- Danby, J.M.A. (1988). *Fundamentals of Celestial Mechanics*. (Standard textbook for orbital mechanics derivations.)
- Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed. (Practical computation methods and lunar elements.)
- Standish, E.M. (1992). "Orbital Ephemerides of the Sun, Moon, and Planets." *Explanatory Supplement to the Astronomical Almanac*. (JPL planetary orbital elements.)
- Montenbruck, O. & Gill, E. (2000). *Satellite Orbits*. (Modern treatment of Kepler's equation and state vector computation.)
