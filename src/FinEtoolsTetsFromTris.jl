module FinEtoolsTetsFromTris

using FinEtools
using TetGen

# tetunsuitable() do pa,pb,pc,pd
#     vol=det(hcat(pb-pa,pc-pa,pd-pa))/6
#     center=0.25*(pa+pb+pc+pd)-[0.5,0.5,0.5]
#     vol> 0.05*norm(center)^2.5
# end

function mesh(fens, bfes)
    input=TetGen.RawTetGenIO{Cdouble}()
    input.pointlist=fens.xyz'
    TetGen.facetlist!(input, connasarray(bfes)')
    output =  tetrahedralize(input, "pQa")
    nfens = deepcopy(fens)
    nfens.xyz = output.pointlist'
    fes = FESetT4(Matrix(FInt.(output.tetrahedronlist')))
    return nfens, fes
end

end # module FinEtoolsTetsFromTris
