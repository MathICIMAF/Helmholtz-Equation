function solve_Helmholtz(f, C)
    
    CPUtic()
    model = GmshDiscreteModel(output_path.*"mesh_geometry.msh") # load geometry and mesh data
    
    Ω = Triangulation(model)

    # Define Lagrangian reference element
    println("Define Lagrangian reference element")
    order = 3;
    reffe = ReferenceFE(lagrangian, Float64, order)
    dΩ = Measure(Ω, 2*order)

    neumanntags = ["Gamma_N"]
    ΓN = BoundaryTriangulation(model,tags=neumanntags)
    dΓN = Measure(ΓN,order)

    sommertags = ["Gamma_R"]
    ΓR = BoundaryTriangulation(model,tags=sommertags)
    dΓR = Measure(ΓR,order)

    # Construct Lagrangian test space with dirichlet boundary condition
    V = TestFESpace(model, reffe, conformity = :H1, dirichlet_tags = ["Gamma_D"], vector_type = Vector{ComplexF64})
    U = TrialFESpace(V, [1]) 
    
    kwav = 2*pi*f/C
    fsource(x) = 0.
    m(u,v) = ∫(kwav^2*(v*u))*dΩ   # mass matrix
    a_bilin(u,v) = ∫( ∇(v)⋅∇(u) - (kwav^2*(v*u)))*dΩ + ∫(1im*kwav*u*v)*dΓR # bilinear form
    b_lin(v) = ∫( v*fsource )*dΩ 
    op = AffineFEOperator(a_bilin,b_lin, U, V)
    beta = 0.5
    SHL(u,v) = 1im*beta*m(u,v)
    opm = AffineFEOperator(SHL,b_lin, U, V)
    A = get_matrix(op)  # Matrix of the system equations
    DF = size(A)
    println("Degree_Freedom= ", DF)
    b = get_vector(op)
    M_shl = get_matrix(opm)

    b = sparse(b)
    A = sparse(A)
    M_shl = sparse(M_shl)
    M_shl = A + M_shl  # Shifter Laplacian Matrix
    println("Computing the ILU Factorization")
    CPUtic()
    fact = ilu(M_shl, τ = 1e-4)
    fact_time = CPUtoc()
    println("fact_time: ", fact_time, " seconds")
    
    println("Computing  the solution")
    DF = DF[1]

    CPUtic()
    u, history = gmres(A, b, Pl = fact, abstol = 1e-5, restart=1, reltol=1e-5, log = true, maxiter = 500) 
    gmres_time = CPUtoc()
    println("GMRES_time: ", gmres_time, " seconds")

    uh = FEFunction(U,u) #building the solution function

    iter = size(history[:resnorm])[1]
    residual = history[:resnorm][end]
    println("min(residual)= ", residual)
    println("iter= ", iter)
    total_time = CPUtoc()
    println("total_time= ", total_time, " seconds")

    data__text ="|\tf(MHz)\t|\tK_Wave\t|\t Degree_Freedom\t|\tIteration\t|\tmin(Residual)\t|\tILU_time(sec)\t|\tGMRES_time(sec)\t|\ttotal_time(sec)\t|\n|\t$f\t|\t$kwav\t|\t$DF\t|\t$iter\t|\t$residual\t|\t$fact_time\t|\t$gmres_time\t|\t$total_time\t|"
    io = open(output_path.*"solucion_data.txt", "w")
    println(io, data__text)
    close(io)

    u_abs = sqrt.(real(u).*real(u)  + imag(u).*imag(u))
    uh_abs = FEFunction(U,complex(u_abs))

    p2 = plot(history[:resnorm], yscale = :log10, label = "GMRES convergency with ILU preconditioning for f=$f MHz", xlabel = "Iteration", ylabel = "Residual norm (preconditioned)", mark = :o)
    display(p2)

    println("saving the solution of Helmholtz equation")
    writevtk(Ω, output_path.*"solucion.vtu", cellfields=["uh_real"=>real(uh),"uh_imag"=>imag(uh), "uh_abs"=>real(uh_abs)])
    
end


println()
println("Solving Helmholtz equation")

solve_Helmholtz(f ,C)
