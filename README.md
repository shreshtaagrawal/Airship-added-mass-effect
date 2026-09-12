# Airship-added-mass-effect
Quantifying how much neglecting added mass changes the predicted longitudinal dynamics of a neutrally buoyant airship.

# Airship Added-Mass Sensitivity Study

## What this is

A comparison of the natural modes of an airship's longitudinal dynamics, computed two ways: once including added mass, and once neglecting it. The point is to measure how much the predictions differ.

## The model

Four states describing motion in the vertical plane:

| State | Symbol | Units |
|---|---|---|
| Altitude | h | m |
| Vertical velocity | v_z | m/s |
| Pitch angle | θ | rad |
| Pitch rate | q | rad/s |

Two kinematic relations, `ḣ = v_z` and `θ̇ = q`, and two dynamic equations from Newton's second law in heave and pitch. These are linearised about level flight and assembled into state-space form,

```
ẋ = A x + B_c δ + B_g w_g
```

where `δ` is the thrust-vector angle and `w_g` is a vertical gust. The natural modes are the eigenvalues of `A`.

## The physics

**Added mass.** Accelerating a body through a fluid means accelerating the surrounding air as well, and that costs force. The entrained mass scales as `k·ρ·V`. Because an airship is neutrally buoyant, its own mass is `m ≈ ρV`, so the two are comparable rather than the added mass being a small correction as it would be for an aircraft. Here it adds 16.5 kg to an 18 kg vehicle.

**Buoyancy pendulum.** The centre of buoyancy sits above the centre of gravity, so tilting the vehicle produces a restoring couple `B·z_b·sin θ`. This acts as a torsional spring and is what makes the pitch mode oscillate. It depends only on geometry, not on airspeed.

**Munk moment.** A slender body in a flow experiences a destabilising pitching moment that grows with the square of airspeed. It opposes the buoyancy pendulum, reducing the net pitch stiffness.

## Results

At a cruise airspeed of 3 m/s:

| Metric | With added mass | Without | Ratio |
|---|---|---|---|
| Pendulum frequency ω_n [rad/s] | 0.7095 | 0.9259 | 1.305 |
| Oscillation period [s] | 9.07 | 7.09 | 0.781 |
| Damping ratio ζ | 0.2181 | 0.2874 | 1.317 |
| Heave time constant τ [s] | 1.438 | 0.766 | 0.533 |

Neglecting added mass makes the vehicle appear 31% faster to oscillate, roughly twice as quick to settle in heave, and 32% better damped. All four errors point the same way: the modelled vehicle is more responsive and better behaved than the real one.

The ratios follow analytically. Rotational quantities scale as `√(I_eff/I_dry) = √(100/60) = 1.29`; the heave time constant scales as `m_dry/m_eff = 18/34.5 = 0.52`. The numerical results differ slightly from these because the full 4×4 system couples pitch and heave, which the isolated formulas ignore.

**Validation.** The model's oscillatory eigenvalue is 0.7095 rad/s against an independently derived pendulum frequency `√(B·z_b/I_eff) = 0.7354 rad/s`, agreeing to within 4%. The small difference is the Munk moment reducing the net stiffness.

## How to run

Open `airship_longitudinal.m` in MATLAB and run it. Base MATLAB only — no toolboxes required. All parameters are in the struct `p` at the top of the script.

## Limitations

- Longitudinal motion only. Lateral and directional dynamics are not modelled.
- Added mass is treated as two scalars, one for heave and one for pitch. The real quantity is a 6×6 tensor, and the off-diagonal coupling terms are ignored.
- The equations are linearised about level flight, so the results hold only for small perturbations. Large-angle behaviour is outside the model's validity.
- The vehicle parameters are representative of a small research airship rather than measured from a specific vehicle.
- The added pitch inertia is estimated rather than derived from the hull geometry.
- Aerodynamic coefficients are approximate values for a streamlined hull with fins.
- The hull is treated as rigid; deformation effects are not included.
