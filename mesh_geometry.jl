println()
println("Creating Geometry and mesh")

function meshsize(elementTag)
    """
    function to calculate the mesh size.
    """
    ElementNode = gmsh.model.mesh.getElementFaceNodes(elementTag[1][1], 4)
    hg = 0
    print("Compute the mesh size")
    for j in 1:4:length(ElementNode)
        n1 = gmsh.model.mesh.getNode(ElementNode[j])
        n2 = gmsh.model.mesh.getNode(ElementNode[j+1])
        n3 = gmsh.model.mesh.getNode(ElementNode[j+2])
        n4 = gmsh.model.mesh.getNode(ElementNode[j+3])

        dist1 = sqrt(sum((n1[1] - n2[1]).^2))
        dist2 = sqrt(sum((n2[1] - n3[1]).^2))
        dist3 = sqrt(sum((n3[1] - n4[1]).^2))
        dist4 = sqrt(sum((n4[1] - n1[1]).^2))           

        hg = maximum([hg, dist1, dist2, dist3, dist4])

    end
    println(", mesh size: ", hg)
    return hg
end

# Generate Mesh and Geometry using Gmsh
gmsh.initialize()

function createGeometryAndMesh()
    CPUtic()
    gmsh.clear()
    gmsh.option.setNumber("General.Terminal", 0);
    gmsh.option.setNumber("Mesh.Algorithm", 8); # (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: FrontalDelaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm",3); # (0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad)
    gmsh.option.setNumber("Geometry.PointNumbers",0);
    gmsh.model.add("semicircle")
    
    geo = gmsh.model.geo;
    
    # Angle to add points to transfinite surface
    α = 45
    
    # Radius of simulation domain
    r = gmsh.onelab.getNumber("Parameters/2Geometry/1Semicircle Radius (r)")[1] 
    # Transducer aperture: 2*a
    a = gmsh.onelab.getNumber("Parameters/2Geometry/2Aperture (a)")[1]  
    println("Radius: ", r," and Aperture: ", a) 
    
    # Transducer aperture: 2*a
    f0 = gmsh.onelab.getNumber("Parameters/1Frequency (f in MHz)")[1]
    f = f0*1e6 # frequency en Hz 
    
    k=2*pi*f/C
    lambda = 2*pi/k    
    
    density = 5;
    mshd_R = r / density; # Mesh density at domain boundary
    mshd_a = a / density; # Mesh density at transducer aperture

    # Points 
    geo.addPoint(0, 0, 0, mshd_a, 1)  # Semicircle boundary midpoint
    geo.addPoint(-r, 0, 0, mshd_R, 2) # Boundary left
    geo.addPoint(-a, 0, 0, mshd_a, 3) # Transducer aperture left
    geo.addPoint(a, 0, 0, mshd_a, 4)  # Transducer aperture right
    geo.addPoint(r, 0, 0, mshd_R, 5)  # Boundary right
    # Additional points for transfinite surface
    geo.addPoint(r*cos(α*pi/180), r*sin(α*pi/180), 0, mshd_R, 8);
    geo.addPoint(r*cos(pi-(α*pi/180)), r*sin(pi-(α*pi/180)), 0, mshd_R, 9);

    # Lines
    geo.addLine(2, 3, 1)
    geo.addLine(3, 4, 2)
    geo.addLine(4, 5, 3)
    # Semicircle
    geo.addCircleArc(5, 1, 8, 4)
    geo.addCircleArc(8, 1, 9, 5)
    geo.addCircleArc(9, 1, 2, 6)

    # Surfaces
    geo.addCurveLoop([1, 2, 3, 4, 5, 6], 1)
    geo.addPlaneSurface([1], 1)

    # Physics Group
    geo.addPhysicalGroup(2, [1], 1)
    geo.addPhysicalGroup(1, [2], 1)
    geo.addPhysicalGroup(1, [1, 3], 2)
    geo.addPhysicalGroup(1, [4,5,6], 3)

    # Region and Boundary Name
    gmsh.model.setPhysicalName(2, 1, "Omega")
    gmsh.model.setPhysicalName(1, 1, "Gamma_D")
    gmsh.model.setPhysicalName(1, 2, "Gamma_N")
    gmsh.model.setPhysicalName(1, 3, "Gamma_R")
   

    geo.synchronize(); # synchronize the view panel

    n = gmsh.onelab.getNumber("Parameters/3Mesh/1PartitionX (n)")[1]
    m = gmsh.onelab.getNumber("Parameters/3Mesh/2PartitionY (m)")[1]
    if n>0 && m>0
        n2 = cld(2*a*n,2*r)
        if n2>n
            n2 = 2
        end      

        if (n+2-n2)%2==0
            n1 = (n+2-n2)/2
        else
            n2 = n2 + 1
            n1 = (n+2-n2)/2
        end

        # Change to Transfinite Curve
        gmsh.model.mesh.setTransfiniteCurve(1,n1);
        gmsh.model.mesh.setTransfiniteCurve(2,n2);
        gmsh.model.mesh.setTransfiniteCurve(3,n1);
        gmsh.model.mesh.setTransfiniteCurve(4,m);
        gmsh.model.mesh.setTransfiniteCurve(5,n);
        gmsh.model.mesh.setTransfiniteCurve(6,m);
        # Make a Transfinite Surface
        gmsh.model.mesh.setTransfiniteSurface(1,"Right",[2,5,8,9]);
        
        
        gmsh.model.mesh.generate(2) # creating the initial mesh
        gmsh.model.mesh.recombine();  # transforming to squares mesh
        gmsh.model.mesh.setRecombine(2,1); # transforming to extructural mesh

        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.renumberNodes()
        gmsh.model.mesh.renumberElements()

        Ref = gmsh.onelab.getNumber("Parameters/3Mesh/3Refine (ref)")[1]
        
        for i in 1:1:Ref
            gmsh.model.mesh.refine()
            gmsh.model.mesh.removeDuplicateNodes()
            gmsh.model.mesh.renumberNodes()
            gmsh.model.mesh.renumberElements()
            n = n*(2)
            m = m*(2)
        end
    
        time_mesh = CPUtoc()
        
        nbcuad = gmsh.option.getNumber("Mesh.NbQuadrangles")
        global nbtriang = gmsh.option.getNumber("Mesh.NbTriangles")
        
        element = gmsh.model.mesh.getElementTypes(2, -1)
        
        hg = meshsize(element)
        
        τ = lambda/hg
        println("λ/h= ", τ)

        println("nbcuad: ", nbcuad,"   nbtriang: ", nbtriang)
        println("time_mesh: ", time_mesh, " sec")
        data__text ="| \t f(MHz) \t | \t K_Wave \t | \t lambda \t | \t r \t | \t a \t | \t n \t |\t m \t |\t h \t | \t λ/h \t |\t time_mesh \t|\n| \t $f \t | \t $k \t | \t $lambda \t | \t $r \t | \t $a \t | \t $n \t | \t $m \t | \t $hg \t | \t $τ \t | \t $time_mesh \t|"
        gmsh.write(output_path.*"mesh_geometry.msh")
        io = open(output_path.*"mesh_data.txt", "w")
        println(io, data__text)
        close(io)
    end
    println()
    println("####################################################")
    println()
    
end

gmsh.onelab.set("""[
    {
        "type":"number",
        "name":"Parameters/1Frequency (f in MHz)",
        "values":[$f0],
        "min":0.25,
        "max":2,
        "step":0.01
    },
    {
        "type":"number",
        "name":"Parameters/2Geometry/1Semicircle Radius (r)",
        "values":[$r],
        "min":0,
        "max":120,
        "step":1
    },
    {
        "type":"number",
        "name":"Parameters/2Geometry/2Aperture (a)",
        "values":[$a],
        "min":0,
        "max":120,
        "step":1
    },
    {
        "type": "number",
        "name": "Parameters/3Mesh/2PartitionY (m)",
        "values":[$mm],
        "min":100,
        "max":1e4,
        "step":50 
    },
    {
        "type": "number",
        "name": "Parameters/3Mesh/1PartitionX (n)",
        "values":[$nn],
        "min":200,
        "max":2e4,
        "step":50 
    },
    {
        "type": "number",
        "name": "Parameters/3Mesh/3Refine (ref)",
        "values":[$ref],
        "min":0,
        "max":5,
        "step":1 
    }    
]""");

# Check for event in view panel
function checkForEvent()
    action = gmsh.onelab.getString("ONELAB/Action")
    if length(action) > 0 && action[1] == "check"
        gmsh.onelab.setString("ONELAB/Action", [""])
        createGeometryAndMesh()
        gmsh.graphics.draw()
    end
    return true
end

createGeometryAndMesh()

if !("-nopopup" in ARGS)
    gmsh.fltk.initialize()
    while gmsh.fltk.isAvailable() == 1 && checkForEvent()
        gmsh.fltk.wait()
    end
    global f0 = gmsh.onelab.getNumber("Parameters/1Frequency (f in MHz)")[1]
    global nn = gmsh.onelab.getNumber("Parameters/3Mesh/1PartitionX (n)")[1]
end

gmsh.finalize();