module FinEtoolsTetsFromTris

using FinEtools
using TetGen
using LinearAlgebra

function _edgestatistics(fens, fes)
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

- `fens`, `bfes`:  finite element node set and a set of the triangles. A check
  is run to see whether the triangles form a closed shape (i.e. the boundary is
  null).
- `tetgen_args`: optional arguments to be passed to `tetrahedralize()`. The
  default arguments controlled the element size by setting the maximum volume.

# Returns

- `nfens`, `fes`: finite element node set and a set of the tetrahedra.
"""
function mesh(fens, bfes; tetgen_args = "")
    input=TetGen.RawTetGenIO{Cdouble}()
    input.pointlist=fens.xyz'
    bbfes = meshboundary(bfes)
    count(bbfes) == 0 || error("The triangular mesh does not appear to be closed")
    TetGen.facetlist!(input, connasarray(bfes)')
    esmin, esmax, esmean = _edgestatistics(fens, bfes)
    maxvol = 1.5 * esmax^3 / 6
    tetgen_args == "" && (tetgen_args = "pQq1.4a$(maxvol)")
    output =  tetrahedralize(input, tetgen_args)
    nfens = FENodeSet(Float64.(output.pointlist'))
    fes = FESetT4(Int.(output.tetrahedronlist'))
    return nfens, fes
end

end # module FinEtoolsTetsFromTris
