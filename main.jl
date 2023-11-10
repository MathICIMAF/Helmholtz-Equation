"""
    HelmholtzEquationCode

    HelmholtzEquation is a FEM code implemented in Julia language to compute the pressure field of a
    Helmholtz problem with mixed boundary conditions.

    ///////////////////////////////////////////////////////////////////////////////////////////
    Input:
    f - frequency (Hz). Use a frequency value in the interval [0.25 , 2]x10^6 Hz.  
    r - Radius of semicircle.
    a - Transducer semiaperture.
    n, m - Mesh partition in X and Y direction.


    Output
    The program save 4 files in the output_file folder. 
    - "mesh_geometry.msh" contains the mesh of the semicircle region.
    - "mesh_data.txt" save in a table the values of someone parameters of the mesh.
    - "solucion.vtu" save the Helmholtz solution to be graph in a view software.
    - "solucion_data.txt" show in a table someone parameters of solver.
"""
# Geometry and Mesh
using Gmsh
using GeometryBasics
# Solve equation
using   Gridap,
        Gridap.Algebra,
        Gridap.FESpaces,
        Gridap.ReferenceFEs,
        Gridap.Arrays,
        Gridap.Geometry,
        Gridap.Fields,
        Gridap.CellData,
        GridapGmsh
using IterativeSolvers
using IncompleteLU
using SparseArrays
using BenchmarkTools 
# graph and time
using Plots
using CPUTime

C = 1500 # Wave velocity in wather
global f0 = 0 #Frequency in MHz

r = 0.133 # Radius of simulation domain
a = 0.01  # Transducer semiaperture: (aperture 2*a)

global nn = 300
mm = 150
ref = 0;

global output_path = "output_file/"
global nbtriang = 0

include.("mesh_geometry.jl")

f = f0*1.e6
println("Frequency: ", f0, " MHz")

if f0<0.25
    println("The frequency must be greater than 0.25 MHz.")
elseif nn<200 && mm<100
    println("Wrong mesh. n must be greater than 200 and m greater than 100.")
elseif nbtriang > 0
    println("Wrong mesh. You have some triangles in your mesh. Try smaller n and m values and use the refine parameter.")
else
    include.("solve_Helmholtz.jl")
end
