# Potential-flow initializer

`potential_flow` builds a Neko field file containing `p`, `u`, `v`, and `w`
from a scalar velocity-potential solve. It is intended to provide a smooth,
approximately divergence-free initial condition before starting an
incompressible `pnpn` simulation.

```console
potential_flow case.case [output.fld]
```

The default output name is `potential_flow.fld`. As with other Neko field
writers, this creates a field series such as `potential_flow0.nek5000` and
`potential_flow0.f00000`. Point the case's field initial condition at the first
data file, for example:

```json
"initial_condition": {
  "type": "field",
  "file_name": "potential_flow0.f00000"
}
```

The utility solves

\[
  \nabla^2 \phi = 0, \qquad \mathbf{u} = \nabla \phi.
\]

The potential equation is solved with GMRES and a PHMG preconditioner. PHMG
uses the settings in `case.fluid.pressure_solver.preconditioner`, while the
convergence tolerance is inherited from
`case.fluid.pressure_solver.absolute_tolerance`. A mesh with no more than one
element per MPI rank falls back to Jacobi because it cannot form a valid PHMG
coarse hierarchy.

Velocity-prescribed boundaries provide the Neumann condition
`d(phi)/dn = u_boundary . n`. Pressure/outflow boundaries set `phi = 0`; the
zero value only fixes the arbitrary potential reference. The discontinuous
element gradients are projected through the assembled mass matrix before the
velocity boundary conditions are reapplied.

Supported boundary mappings are:

- Neumann potential plus velocity reapplication: `velocity_value`,
  `expression_velocity`, `blasius_profile`, stationary `no_slip`, and
  `symmetry`.
- Homogeneous Neumann potential: `shear_stress` and `wall_model`.
- Zero Dirichlet potential: pressure and outflow variants.

Moving walls, `user_velocity`, and `overset_interface` require runtime user
data and are rejected. CPU and configured Neko device backends are supported;
the potential solve, gradient projection, boundary application, and diagnostics
remain on the device until the result is copied to the host for output. For a
pure-Neumann problem, the prescribed total flux must be compatible; the utility
diagnoses and rejects a nonzero net flux. A fully periodic domain has no way to
infer a nonzero mean velocity and therefore produces the zero-flow solution.

The reported integrated prescribed Neumann flux only includes boundaries where
the normal velocity is prescribed. It need not be zero when the problem also
has Dirichlet-potential boundaries: their normal flux is part of the solution
and balances the prescribed contribution.

The field pressure is initialized to zero. The utility reports the potential
solver residual, integrated prescribed Neumann flux, and velocity-divergence
norms before and after strong velocity boundary conditions are reapplied.
