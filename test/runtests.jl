module m001
using Test
using FinEtools
using FinEtools.MeshExportModule.VTKWrite
using FinEtoolsTetsFromTris

function test()
    fens, fes = T4block(1.0, 2.0, 3.0, 11, 22, 33)
    bfes = meshboundary(fes)
    connected = findunconnnodes(fens, bfes);
    fens, new_numbering = compactnodes(fens, connected);
    bfes = renumberconn!(bfes, new_numbering);
    fens, fes = FinEtoolsTetsFromTris.mesh(fens, bfes)
    VTKWrite.vtkwrite("m001.vtu", fens, fes)
    true
end

test()
nothing
end
