module FinEtoolsTetsFromTris

using FinEtools
using TetGen
using LinearAlgebra


"""
    edgestatistics(fens, fes)

Calculate statistics of the edges of the surface mesh.
"""
function edgestatistics(fens, fes)
    ncnodes = nodesperelem(fes)
    esmin = Inf
    esmax = 0.0
    esmean = 0.0
    nedges = 0
    for i in 1:count(fes)
        for f in 1:ncnodes
            for s in f+1:ncnodes
                es = norm(fens.xyz[fes.conn[i][s], :] - fens.xyz[fes.conn[i][f], :])
                esmin = min(esmin, es)
                esmax = max(esmax, es)
                esmean += es
                nedges += 1
            end
        end
    end
    return esmin, esmax, esmean / nedges
end

"""
    mesh(fens, bfes; tetgen_args = "")

Generate tetrahedra to fill volume bounded by triangles.

# Arguments

- `fens`, `bfes`:  finite element node set and a set of the boundary triangles.
  A check is run to see whether the triangles form a closed shape (i.e. the
  boundary of the boundary is null).
- `tetgen_args`: optional arguments to be passed to `tetrahedralize()`. The
  default argument makes the tetrahedron generation quiet (`"Q"`).

# Returns

- `nfens`, `fes`: finite element node set and a set of the tetrahedra.

# Notes

Constraints can be added to tetrahedralization by passing additional
arguments to `tetrahedralize()`. For example, to limit the maximum volume of 
the tetrahedra, you can use the following code:   
```
    esmin, esmax, esmean = edgestatistics(fens, bfes)
    maxvol = 1.5 * esmax^3 / 6
    nfens, fes = mesh(fens, bfes; tetgen_args = "Qa\$(maxvol)")
```
"""
function mesh(fens, bfes; tetgen_args = "pQq1.4")
    input=TetGen.RawTetGenIO{Cdouble}()
    input.pointlist=fens.xyz'
    bbfes = meshboundary(bfes)
    count(bbfes) == 0 || error("The triangular mesh does not appear to be closed")
    TetGen.facetlist!(input, connasarray(bfes)')
    output =  tetrahedralize(input, tetgen_args)
    nfens = FENodeSet(Float64.(output.pointlist'))
    fes = FESetT4(Int.(output.tetrahedronlist'))
    return nfens, fes
end

end # module FinEtoolsTetsFromTris
