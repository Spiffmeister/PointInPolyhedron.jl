# PointInPolygon.jl


Julia implementation of two methods for detecting if a point lies inside a triangulated volume.




```julia
] add https://github.com/Spiffmeister/PointInPolygon.git
```

```julia
using PointInPolygon
```


```julia
include("test/torus.jl")
```



## A quick example

As an example consider the unit cube

```julia
unit_cube_vertex = Tuple(SVector(i,j,k) for i ∈ [0.0,1.0], j ∈ [0.0,1.0], k ∈ [0.0,1.0])
connections = ()

m = Mesh{Float64}(vertices, connections)
```

Then given a set of points

```julia
npts = 10

pts = Tuple(randn(3) for _ in 1:npts)
```

The winding number can be computed

Will return `0` if the point is inside the domain and `1` if outside.

```julia
winding_number(m, pts)
```
