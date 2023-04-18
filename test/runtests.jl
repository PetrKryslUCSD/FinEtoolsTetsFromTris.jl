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

module m002
using Test
using FinEtools
using FinEtools.MeshExportModule.VTKWrite
using FinEtools.MeshImportModule
using FinEtoolsTetsFromTris

function test()
    output = MeshImportModule.import_ABAQUS(
        dirname(@__FILE__) * "/" * "spherical_shell-0_030.inp";
    )
    fens = output["fens"]
    fes = nothing
    for j in 1:length(output["fesets"])
        if output["fesets"][j] isa FESetT4
            fes1 = output["fesets"][j]
            setlabel!(fes1, j)
            if fes == nothing
                fes = fes1
            else
                fes = cat(fes, fes1)
            end
        end
    end

    # @test count(output["fens"]) == 1406
    # @test count(output["fesets"][1]) == 829

    bfes = meshboundary(fes)
    connected = findunconnnodes(fens, bfes);
    fens, new_numbering = compactnodes(fens, connected);
    bfes = renumberconn!(bfes, new_numbering);
    fens, fes = FinEtoolsTetsFromTris.mesh(fens, bfes)
    VTKWrite.vtkwrite("m002.vtu", fens, fes)
    true
end

test()
nothing
end

module m003
using Test
using FinEtools
using FinEtools.MeshExportModule.VTKWrite
using FinEtools.MeshImportModule
using FinEtoolsTetsFromTris
using LinearAlgebra

function meshstatistics(fens, fes)
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

function test()
    output = MeshImportModule.import_ABAQUS(
        dirname(@__FILE__) * "/" * "spherical_shell-0_030.inp";
    )
    fens = output["fens"]
    fes = nothing
    for j in 1:length(output["fesets"])
        if output["fesets"][j] isa FESetT4
            fes1 = output["fesets"][j]
            setlabel!(fes1, j)
            if fes == nothing
                fes = fes1
            else
                fes = cat(fes, fes1)
            end
        end
    end

    # @test count(output["fens"]) == 1406
    # @test count(output["fesets"][1]) == 829

    bfes = meshboundary(fes)
    connected = findunconnnodes(fens, bfes);
    fens, new_numbering = compactnodes(fens, connected);
    bfes = renumberconn!(bfes, new_numbering);
    esmin, esmax, esmean = meshstatistics(fens, bfes)
    maxvol = 1.5 * esmax^3 / 6 * 10
    fens, fes = FinEtoolsTetsFromTris.mesh(fens, bfes; tetgen_args = "pQq1.4a$(maxvol)")
    VTKWrite.vtkwrite("m003.vtu", fens, fes)
    true
end

test()
nothing
end
