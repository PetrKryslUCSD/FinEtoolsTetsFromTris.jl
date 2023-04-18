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

function mesh(fens, bfes; tetgen_args = "")
    input=TetGen.RawTetGenIO{Cdouble}()
    input.pointlist=fens.xyz'
    TetGen.facetlist!(input, connasarray(bfes)')
    esmin, esmax, esmean = _edgestatistics(fens, bfes)
    maxvol = 1.5 * esmax^3 / 6
    tetgen_args == "" && (tetgen_args = "pQq1.4a$(maxvol)")
    output =  tetrahedralize(input, tetgen_args)
    nfens = deepcopy(fens)
    nfens.xyz = output.pointlist'
    fes = FESetT4(output.tetrahedronlist')
    return nfens, fes
end

end # module FinEtoolsTetsFromTris
