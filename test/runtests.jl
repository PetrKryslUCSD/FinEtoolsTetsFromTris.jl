module m00x01
using Test
using FinEtools
using FinEtools.MeshExportModule.VTKWrite
using FinEtools.MeshExportModule.NASTRAN
using FinEtools.MeshImportModule
using FinEtoolsTetsFromTris
using DelimitedFiles

function test()
    xyz = readdlm("uxoxyz.txt", ',')
    fens = FENodeSet(xyz)
    conn = readdlm("uxoconn.txt", ',')
    bfes = FESetT3(Int.(conn))
    fens, fes = FinEtoolsTetsFromTris.mesh(fens, bfes)

    VTKWrite.vtkwrite("m00x01.vtu", fens, fes)

    e = NASTRANExporter("uxo.nas")
    BEGIN_BULK(e)
    for j in eachindex(fens)
        GRID(e, j, fens.xyz[j, :])
    end
    PSOLID(e, 1, 1)
    conn = connasarray(fes)
    for j in eachindex(fes)
        CTETRA(e, j, 1, conn[j, :])
    end
    ENDDATA(e)
    CEND(e)
    close(e)

    output = MeshImportModule.import_NASTRAN(
        "uxo.nas"
    )
    @test count(output["fens"]) == count(fens)
    @test count(output["fesets"][1]) == count(fes)
    true
end

test()
nothing
end
module m00x02
using Test
using FinEtools
using FinEtools.MeshExportModule.VTKWrite
using FinEtools.MeshExportModule.NASTRAN
using FinEtools.MeshImportModule
using FinEtoolsTetsFromTris
using DelimitedFiles

function test()
    xyz = readdlm("uxoxyz.txt", ',')
    fens = FENodeSet(xyz)
    conn = readdlm("uxoconn.txt", ',')
    bfes = FESetT3(Int.(conn))
    esmin, esmax, esmean = FinEtoolsTetsFromTris.edgestatistics(fens, bfes)
    maxvol = 1.5 * esmax^3 / 6
    fens, fes = FinEtoolsTetsFromTris.mesh(fens, bfes; tetgen_args = "Qpq1.4a$(maxvol)")
    
    VTKWrite.vtkwrite("m00x02.vtu", fens, fes)

    e = NASTRANExporter("uxo.nas")
    BEGIN_BULK(e)
    for j in eachindex(fens)
        GRID(e, j, fens.xyz[j, :])
    end
    PSOLID(e, 1, 1)
    conn = connasarray(fes)
    for j in eachindex(fes)
        CTETRA(e, j, 1, conn[j, :])
    end
    ENDDATA(e)
    CEND(e)
    close(e)

    output = MeshImportModule.import_NASTRAN(
        "uxo.nas"
    )
    @test count(output["fens"]) == count(fens)
    @test count(output["fesets"][1]) == count(fes)
    true
end

test()
nothing
end

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
    @test count(fes) == 42728
    VTKWrite.vtkwrite("m002.vtu", fens, fes)
    true
end

test()
nothing
end

module m002b
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
    
    esmin, esmax, esmean = FinEtoolsTetsFromTris.edgestatistics(fens, bfes)
    maxvol = 1.5 * esmax^3 / 6
    fens, fes = FinEtoolsTetsFromTris.mesh(fens, bfes; tetgen_args = "Qpq1.4a$(maxvol)")
    
    @test count(fes) == 64269
    VTKWrite.vtkwrite("m002b.vtu", fens, fes)
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
    @test count(fes) == 42928
    VTKWrite.vtkwrite("m003.vtu", fens, fes)
    true
end

test()
nothing
end

module m004
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
    # create an invalid surface set
    bfes = subset(bfes, [1, 3])
    @test_throws ErrorException fens, fes = FinEtoolsTetsFromTris.mesh(fens, bfes)

    true
end

test()
nothing
end



module m005
using Test
using FinEtools
using FinEtools.MeshExportModule.VTKWrite
using FinEtoolsTetsFromTris

function test()
    fens, fes = T4quartercyln(1.0, 2.0, 7, 9; orientation = :b)
    bfes = meshboundary(fes)
    connected = findunconnnodes(fens, bfes);
    fens, new_numbering = compactnodes(fens, connected);
    bfes = renumberconn!(bfes, new_numbering);
    fens, fes = FinEtoolsTetsFromTris.mesh(fens, bfes)
    VTKWrite.vtkwrite("m005.vtu", fens, fes)
    true
end

test()
nothing
end
