using StaticArrays
using PointInPolyhedron
using Tests

# Build the unit cube
unit_cube_vertex = Tuple(SVector(i, j, k) for i ∈ [0.0, 1.0], j ∈ [0.0, 1.0], k ∈ [0.0, 1.0])
connections = ((1, 2, 3), (2, 3, 4))
m = Mesh{Float64}(unit_cube_vertex, connections)

# Generate the test points
test_points = Tuple(randn(3) for _ in 1:100)

# Detect if inside the unit cube
point_inside = [all(0.0 .≤ pt .≤ 1.0) for pt in test_points]

# Solid angle check
sa = solid_angle(m, test_points)
@test all(sa .== point_inside)

# Winding number check
wn = winding_number(m, test_points)
@test all(wn .== point_inside)
