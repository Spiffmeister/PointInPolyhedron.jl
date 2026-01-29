# PointInPolygon.jl


Julia implementation of two methods for detecting if a point lies inside a triangulated volume.




```julia
] add https://github.com/Spiffmeister/PointInPolyhedron.git
```




```julia
include("test/torus.jl")
```



## A quick example



As an example consider the unit cube

```julia
using PointInPolyhedron
using StaticArrays

unit_cube_vertex = Tuple(SVector(i,j,k) for i ∈ [0.0,1.0], j ∈ [0.0,1.0], k ∈ [0.0,1.0])
connections = (
    (2, 1, 3), (3, 4, 2), (6, 2, 4), (4, 8, 6),
    (5, 6, 8), (8, 7, 5), (1, 5, 7), (7, 3, 1),
    (4, 3, 7), (7, 8, 4), (6, 5, 1), (1, 2, 6),
)

m = Mesh{Float64}(vertices, connections)
```

Then given a set of points

```julia
npts = 10

pts = Tuple(randn(3) for _ in 1:npts)
```

We can test if they are inside or outside the volume either by computing the solid angle,
```julia
solid_angle(m, pts)
```
which returns zero if the point it outside the domain and one otherwise.

Alternatively by using the winding number,
```julia
winding_number(m, pts)
```
which returns the same.
