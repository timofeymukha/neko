# Potential-flow initializer

`potential_flow` builds a Neko field file containing `p`, `u`, `v`, and `w`
from a scalar velocity-potential solve. It is intended to provide a smooth,
approximately divergence-free initial condition before starting an
incompressible `pnpn` simulation.

```console
potential_flow potential_flow.json
```

The argument is a standalone initializer configuration. It references the Neko
case but owns the potential solver settings, so changing the simulation's
pressure solver does not change the generated initial condition.

```json
{
  "version": 1.0,
  "case_file": "pseudo2d_5e.json",
  "mesh_file": "mesh.nmsh",
  "polynomial_order": 3,
  "output_filename": "potential_flow.fld",
  "evaluation_time": 0.0,
  "solver": {
    "type": "gmres",
    "preconditioner": {
      "type": "phmg"
    },
    "absolute_tolerance": 1e-6,
    "max_iterations": 500,
    "monitor": true
  }
}
```

`version`, `case_file`, and `solver` are required. The `solver` object uses the
same schema as a Neko linear solver. In particular, all PHMG or HSMG parameters
are placed inside `solver.preconditioner` in the usual form. If omitted,
`solver.max_iterations` defaults to 5000 and `solver.monitor` to `false`.
Absolute and optional relative tolerances are both honored. Projection-space
settings are not used because the initializer performs only one solve.

The following values can be omitted and are then inherited from the referenced
case: `mesh_file`, `polynomial_order`, `cyclic`, and `boundary_conditions`.
`evaluation_time` defaults to the case start time and controls the evaluation
of time-dependent boundary expressions. `output_filename` defaults to
`potential_flow.fld`. Paths written in the initializer configuration are
relative to that file; paths inherited from the Neko case remain relative to
the case file.

An optional `boundary_conditions` array replaces the complete
`case.fluid.boundary_conditions` array. It uses the standard Neko pnpn boundary
condition schema, though the unsupported runtime-user boundary types listed
below are still rejected by the utility.

The configuration can be checked against
`doc/schemas/potential-flow.schema.json`:

```console
python contrib/validate_case_schema.py \
  --schema doc/schemas/potential-flow.schema.json potential_flow.json
```

As with other Neko field writers, an output name of `potential_flow.fld`
creates a field series such as `potential_flow0.nek5000` and
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

The potential equation uses the solver and preconditioner selected in the
initializer configuration. The available scalar Krylov solvers and
preconditioners are the same backend-specific implementations used by Neko.

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
