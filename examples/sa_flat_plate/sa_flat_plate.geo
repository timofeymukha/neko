//
// 3D flat-plate mesh for SA verification (TMR-style domain)
// Domain: x in [-1, 2], y in [0, 1], plate on y=0 for x in [0, 2]
//
// This is an SEM-oriented element mesh (not point-for-point FD/FV grid).
// Chosen as a medium-size SEM mesh: 34 x 24 elements.
// A practical mapping is with lx=5 (4 GLL intervals/element), which gives
// 137 x 97 GLL points overall.
//

SetFactory("OpenCASCADE");

// -------------------------
// Geometry parameters
// -------------------------
x_min = DefineNumber[ -1.0, Name "Domain/x_min" ];
x_le  = DefineNumber[  0.0, Name "Domain/x_le"  ];
x_max = DefineNumber[  2.0, Name "Domain/x_max" ];
y_min = DefineNumber[  0.0, Name "Domain/y_min" ];
y_max = DefineNumber[  1.0, Name "Domain/y_max" ];
z_len = DefineNumber[  0.10, Name "Domain/z_len" ];

// -------------------------
// Element counts (medium)
// -------------------------
// Upstream of LE:  6 elements
// Plate section:  28 elements
// Wall-normal:    24 elements
nx_up    = DefineNumber[  6, Name "Mesh/nx_up"    ];
nx_plate = DefineNumber[ 28, Name "Mesh/nx_plate" ];
ny       = DefineNumber[ 24, Name "Mesh/ny"       ];

// -------------------------
// Stretching controls
// -------------------------
// Cluster near leading edge from both sides.
rx_up    = DefineNumber[ 1.20, Name "Stretch/rx_up"    ];
rx_plate = DefineNumber[ 1.06, Name "Stretch/rx_plate" ];

// Cluster toward wall (y=0).
ry_wall  = DefineNumber[ 1.18, Name "Stretch/ry_wall"  ];

// -------------------------
// Points
// -------------------------
p1 = newp; Point(p1) = {x_min, y_min, 0.0, 1.0};
p2 = newp; Point(p2) = {x_le , y_min, 0.0, 1.0};
p3 = newp; Point(p3) = {x_max, y_min, 0.0, 1.0};
p4 = newp; Point(p4) = {x_min, y_max, 0.0, 1.0};
p5 = newp; Point(p5) = {x_le , y_max, 0.0, 1.0};
p6 = newp; Point(p6) = {x_max, y_max, 0.0, 1.0};

// -------------------------
// Curves
// -------------------------
// Bottom boundary split at leading edge.
l_bottom_up    = newl; Line(l_bottom_up)    = {p1, p2};
l_bottom_plate = newl; Line(l_bottom_plate) = {p2, p3};

// Top boundary split to match block structure.
l_top_up       = newl; Line(l_top_up)       = {p4, p5};
l_top_plate    = newl; Line(l_top_plate)    = {p5, p6};

// Vertical boundaries.
l_left         = newl; Line(l_left)         = {p1, p4};
l_le           = newl; Line(l_le)           = {p2, p5};
l_right        = newl; Line(l_right)        = {p3, p6};

// -------------------------
// Surfaces (two blocks)
// -------------------------
ll_up = newll;
Curve Loop(ll_up) = {l_bottom_up, l_le, -l_top_up, -l_left};
s_up = news;
Plane Surface(s_up) = {ll_up};

ll_plate = newll;
Curve Loop(ll_plate) = {l_bottom_plate, l_right, -l_top_plate, -l_le};
s_plate = news;
Plane Surface(s_plate) = {ll_plate};

// -------------------------
// Structured meshing
// -------------------------
n_up_pts = nx_up + 1;
n_plate_pts = nx_plate + 1;
n_y_pts = ny + 1;

// x-clustering near the leading edge (x=0) from both sides.
Transfinite Curve {l_bottom_up, l_top_up} = n_up_pts Using Progression
  (1.0 / rx_up);
Transfinite Curve {l_bottom_plate, l_top_plate} = n_plate_pts
  Using Progression rx_plate;

// y-clustering toward wall (y=0) for all vertical lines.
Transfinite Curve {l_left, l_le, l_right} = n_y_pts Using Progression
  ry_wall;

Transfinite Surface {s_up};
Transfinite Surface {s_plate};
Recombine Surface {s_up, s_plate};

// -------------------------
// 3D extrusion: one element in z
// -------------------------
out_up[] = Extrude {0.0, 0.0, z_len} {
  Surface{s_up};
  Layers{1};
  Recombine;
};

out_plate[] = Extrude {0.0, 0.0, z_len} {
  Surface{s_plate};
  Layers{1};
  Recombine;
};

// Volumes
Physical Volume("fluid_volume", 10) = {out_up[1], out_plate[1]};

// 3D boundaries
Physical Surface("inlet", 2) = {out_up[3]};
Physical Surface("outlet", 3) = {out_plate[3]};
Physical Surface("farfield", 4) = {out_up[4], out_plate[4]};
Physical Surface("wall", 5) = {out_plate[2]};
Physical Surface("symmetry", 6) = {out_up[2]};

// z-normal planes (front/back)
Physical Surface("z_minus", 7) = {s_up, s_plate};
Physical Surface("z_plus", 8) = {out_up[0], out_plate[0]};

Coherence;

// -------------------------
// Mesh settings
// -------------------------
Mesh.MshFileVersion = 2.2;
Mesh.RecombineAll = 1;
Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll = 0;
Mesh.Binary = 0;

SetOrder 2;

Save "sa_flat_plate.msh";
