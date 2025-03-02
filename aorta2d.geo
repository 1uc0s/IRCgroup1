// Parameters for aorta with atherosclerosis
diameter = 25.0;  // Aorta diameter in mm (average adult aorta ~25mm)
length = 150.0;   // Length of segment in mm
stenosis_level = 0.5;  // 0.0 = no stenosis, 1.0 = complete occlusion (adjust as needed)
stenosis_length = 20.0;  // Length of the stenotic region in mm
stenosis_position = length/2;  // Position of stenosis center from inlet
mesh_size = 0.5;  // Default mesh size (smaller values create finer mesh)

// Calculate stenosis height based on stenosis level
stenosis_height = diameter * stenosis_level / 2;

// Points for upper wall
Point(1) = {0, diameter/2, 0, mesh_size};  // Inlet top
Point(2) = {stenosis_position - stenosis_length/2, diameter/2, 0, mesh_size};  // Start of stenosis top
Point(3) = {stenosis_position, diameter/2 - stenosis_height, 0, mesh_size/4};  // Peak of stenosis top (finer mesh)
Point(4) = {stenosis_position + stenosis_length/2, diameter/2, 0, mesh_size};  // End of stenosis top
Point(5) = {length, diameter/2, 0, mesh_size};  // Outlet top

// Points for lower wall
Point(6) = {0, -diameter/2, 0, mesh_size};  // Inlet bottom
Point(7) = {stenosis_position - stenosis_length/2, -diameter/2, 0, mesh_size};  // Start of stenosis bottom
Point(8) = {stenosis_position, -diameter/2 + stenosis_height, 0, mesh_size/4};  // Peak of stenosis bottom (finer mesh)
Point(9) = {stenosis_position + stenosis_length/2, -diameter/2, 0, mesh_size};  // End of stenosis bottom
Point(10) = {length, -diameter/2, 0, mesh_size};  // Outlet bottom

// Create splines for the walls (smoother than straight lines)
Spline(1) = {1, 2, 3, 4, 5};  // Upper wall
Spline(2) = {6, 7, 8, 9, 10};  // Lower wall

// Create inlet and outlet lines
Line(3) = {1, 6};  // Inlet
Line(4) = {5, 10}; // Outlet

// Create surface for the fluid domain
Line Loop(1) = {1, 4, -2, -3};
Plane Surface(1) = {1};

// Define physical groups for OpenFOAM boundary conditions
Physical Curve("inlet") = {3};
Physical Curve("outlet") = {4};
Physical Curve("wall") = {1, 2};
Physical Surface("fluid") = {1};

// Mesh settings for better quality
Mesh.RecombineAll = 1;  // Generate quadrilateral elements where possible
Mesh.Smoothing = 20;    // Mesh smoothing steps
Mesh.Algorithm = 6;     // Frontal-Delaunay for quads