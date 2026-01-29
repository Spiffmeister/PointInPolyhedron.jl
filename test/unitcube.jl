using StaticArrays
using PointInPolyhedron
using Random
using Test

# Build the unit cube
unit_cube_vertex = Tuple(SVector(i, j, k) for i ∈ [0.0, 1.0], j ∈ [0.0, 1.0], k ∈ [0.0, 1.0])

# cube = Polyhedron(
#     vertex_positions=[
#         (-1, -1, -1), (-1, -1, +1), (-1, +1, -1), (-1, +1, +1),
#         (+1, -1, -1), (+1, -1, +1), (+1, +1, -1), (+1, +1, +1),
#     ],
#     triangles=[
#         [1, 3, 2], [1, 0, 4], [1, 5, 7],
#         [2, 0, 1], [2, 6, 4], [2, 3, 7],
#         [4, 5, 1], [4, 0, 2], [4, 6, 7],
#         [7, 3, 1], [7, 6, 2], [7, 5, 4],
#     ],
# )

connections = (
    # 000 - 100 - 110 - 010 -- Front
    (2, 1, 3),
    (3, 4, 2),
    # 100 - 101 - 111 - 110 -- Right
    (6, 2, 4),
    (4, 8, 6),
    # 101 - 111 - 011 - 001 -- Back
    (5, 6, 8),
    (8, 7, 5),
    # 001 - 011 - 010 - 000 -- Left
    (1, 5, 7),
    (7, 3, 1),
    # 010 - 110 - 111 - 011 -- Top
    (4, 3, 7),
    (7, 8, 4),
    # 001 - 000 - 100 - 101 -- Bottom
    (6, 5, 1),
    (1, 2, 6),
)

m = Mesh{Float64}(unit_cube_vertex, connections)

# Generate the test points
Random.seed!(123)
test_points = Tuple(rand(3) * 1.5 for _ in 1:100)

# Detect if inside the unit cube
point_inside = Tuple(Float64(all(0.0 .≤ pt .≤ 1.0)) for pt in test_points)

# Solid angle check
sa = solid_angle(m, test_points)
sa_inout = Float64.(sa .> 0.5)
@test all(sa_inout .== point_inside)

# Winding number check
wn = winding_number(m, test_points)
@test all(wn .== point_inside)
