using Molly
using Printf
using EzXML
include("AtomParameters.jl")
include("generateFragment.jl")

#This file contains utilities for generating structures for use with Molly, as well as structural analysis

#This is a string example of an .xyz file as produced by Orca or Chemcraft or Avogadro
#makeAtoms should be able to read this string and produce a list of atoms
# fileString = "O     0.000000000000      0.000000000000      0.000000000000
#               H    -0.000000000008      0.759336999999      0.596043000001
#               H    -0.000000000008     -0.759337000001      0.596042999999"

function makeAtoms(inputXYZ)

    lines = split(inputXYZ,"\n")
    
    if isa(tryparse(Int,lines[1]),Int) #if the first line is a number, it means that it is an xmol file (there are two lines of header)
        n_atoms = parse(Int,lines[1])
        shift_k = 2
    else 
        n_atoms = size(lines)[1]     #the file just starts with the atom labels and coordinates, so I just count the lines
        shift_k = 0
    end
    #Maybe here there is a problem with using type Any for Atoms. 
    #atoms = Vector{Any}(undef, n_atoms)
    #atoms = Vector{Molly.Atom}(undef)
    #atoms = fill(Atom(mass=0.0),0)
    atoms = []
    #This is a clean coordinate matrix, generally in Angstrom (but could be Bohr)
    coordsXYZ = Matrix{Float64}(undef, n_atoms, 3)
    

    for k = (1:n_atoms).+shift_k
        line_contents = split(lines[k])
        #I need to make a dictionary for Lennard-Jones parameters of atoms
        elementType = line_contents[1]
        atomProperties = elementTypeDict[elementType]
        atomMass = atomProperties[2]
        atomLJ_σ = atomProperties[3]
        atomLJ_ϵ = atomProperties[4]
        push!(atoms,Molly.Atom(mass=(atomMass)u"g/mol", σ=(atomLJ_σ)u"nm", ϵ=(atomLJ_ϵ)u"kJ * mol^-1"))

        coordsXYZ[k,:] = parse.(Float64,line_contents[2:4])./10 #to convert from angstrom to nm
    end

    #print(coordsXYZ)
    #coords = []
    #println("Initialized coords")
    #println(coords)
    #This is not really used, becuase place_atoms_at_coordinates multiplies the boundary by zero. It's a hack
    boundary = CubicBoundary((3.0)u"nm")
    coords = place_atoms_at_coordinates(n_atoms, coordsXYZ, boundary)

    #This is still not working properly, need to use place_atoms_at_coordinates 
    # for k = 1:n_atoms
    #     #println("printing SVector of coordinates")
    #     #println(SVector{3}(coordsXYZ[k,1],coordsXYZ[k,2],coordsXYZ[k,3])u"nm")
    #     push!(coords, [(coordsXYZ[k,1],coordsXYZ[k,2],coordsXYZ[k,3])u"nm"])
    #     #push!(coords, SVector{3}(coordsXYZ[k,1],coordsXYZ[k,2],coordsXYZ[k,3])u"nm")

    # end

    return atoms, coords, coordsXYZ

end

        function random_coord(boundary::CubicBoundary{3, T}; rng=Random.default_rng()) where T
            return rand(rng, SVector{3, T}) .* boundary
        end


#This function is (not) the same (anymore) as place_atoms, but it uses XYZ coordinates from an array
function place_atoms_at_coordinates(n_atoms::Integer,
                     coordsXYZ,
                     boundary)
    rng=Random.default_rng()
    coords = SArray[]
    sizehint!(coords, n_atoms)

    ####This is what changes from place_atoms
    iXYZ = 1
    ####
    while length(coords) < n_atoms
        ####This is what changes from place_atoms
        new_random_coord = random_coord(boundary; rng=rng).*0
        # print("iXYZ: ")
        # println(iXYZ)
        # println("the new_random_coord is")
        # println(new_random_coord)
        new_coord = new_random_coord.*0 + (coordsXYZ[iXYZ,:]-[ustrip(boundary)[1], ustrip(boundary)[2], ustrip(boundary)[3]]./2)u"nm"
        # println("the new_coord is")
        # println(new_coord)
        iXYZ = iXYZ + 1
        ####
        okay = true

        if okay
            push!(coords, new_coord)
        end
    end
    return [coords...]
end

#3/09/2026: In the makeBonds, makeAngles, makeDihedrals Functions
#I should add as input the atomTypes, as well as the forcefield dictionaries
#in order to find and get the parameters from the ff

function makeBonds(atoms,topology,atomTypeArray,BondTypes)
    #for bonds, we only need topology.bondIndList
    n_atoms = size(atoms,1)
    n_bonds_tot = size(topology.bondIndList,1)
    
    Interaction = []
    for j in 1:n_bonds_tot
        id1 = Int64(topology.bondIndList[j,1])
        id2 = Int64(topology.bondIndList[j,2])

        element1 = atomTypeArray[id1]["element"]
        hyb1 = atomTypeArray[id1]["hybridization"]
        class1 = atomTypeArray[id1]["class"]
        gafftype1 = atomTypeArray[id1]["gafftype"]

        element2 = atomTypeArray[id2]["element"]
        hyb2 = atomTypeArray[id2]["hybridization"]
        class2 = atomTypeArray[id2]["class"]
        gafftype2 = atomTypeArray[id2]["gafftype"]

        selBondTypes = filter(d -> (d["type1"] == gafftype1 && d["type2"] == gafftype2) || (d["type1"] == gafftype2 && d["type2"] == gafftype1), BondTypes)
        # println(id1)
        # println(gafftype1)
        # println(id2)
        # println(gafftype2)
        # println(selBondTypes)

        #k = 200_000u"kJ * mol^-1 * nm^-2"

        # if (element1 == "C")
        #     if (element2 == "C")
        #         k = 300_000u"kJ * mol^-1 * nm^-2"
        #     elseif (element2 == "O")
        #         k = 350_000u"kJ * mol^-1 * nm^-2"
        #     end
        # end

        # if (hyb1 == "sp2") && (hyb2 == "sp2")
        #     k = k*1.5
        # end
        if !isempty(selBondTypes)
            #println(selBondTypes)
            #println(selBondTypes[1]["k"])
            push!(Interaction,HarmonicBond(k=(selBondTypes[1]["k"])u"kJ * mol^-1 * nm^-2", r0=(selBondTypes[1]["length"])u"nm"))
        else
            error("Bond parameters for atoms $gafftype1 and $gafftype2 were not found in BondTypes.")
        end
    end


    bonds = InteractionList2Atoms(
    Vector{Int64}(topology.bondIndList[:,1]),           # First atom indices
    Vector{Int64}(topology.bondIndList[:,2]),           # Second atom indices
    Interaction
    )
    return bonds
end

function makeAngles(atoms,topology,atomTypeArray,AngleTypes)
    n_atoms = size(atoms,1)
    n_angles_tot = size(topology.angleIndList,1)
    
    Interaction = []
    for j in 1:n_angles_tot
        id1 = Int64(topology.angleIndList[j,1])
        id2 = Int64(topology.angleIndList[j,2])
        id3 = Int64(topology.angleIndList[j,3])

        element1 = atomTypeArray[id1]["element"]
        hyb1 = atomTypeArray[id1]["hybridization"]
        class1 = atomTypeArray[id1]["class"]
        gafftype1 = atomTypeArray[id1]["gafftype"]

        element2 = atomTypeArray[id2]["element"]
        hyb2 = atomTypeArray[id2]["hybridization"]
        class2 = atomTypeArray[id2]["class"]
        gafftype2 = atomTypeArray[id2]["gafftype"]

        element3 = atomTypeArray[id3]["element"]
        hyb3 = atomTypeArray[id3]["hybridization"]
        class3 = atomTypeArray[id3]["class"]
        gafftype3 = atomTypeArray[id3]["gafftype"]

        #For angles we keep the central atom fixed (type2), and allow both possibilities for type1 and type3
        selAngleTypes = filter(d -> (d["type1"] == gafftype1 && d["type2"] == gafftype2 && d["type3"] == gafftype3) || (d["type1"] == gafftype3 && d["type2"] == gafftype2 && d["type3"] == gafftype1), AngleTypes)
        # println(id1)
        # println(gafftype1)
        # println(id2)
        # println(gafftype2)
        # println(selBondTypes)

        if !isempty(selAngleTypes)
            #println(selBondTypes)
            #println(selBondTypes[1]["k"])
            push!(Interaction,HarmonicAngle(k=(selAngleTypes[1]["k"])u"kJ/mol", θ0=(selAngleTypes[1]["angle"])u"rad"))
        else
            error("Angle parameters for atoms $gafftype1, $gafftype2 and $gafftype3 were not found in AngleTypes.")
        end
    end

    angles = InteractionList3Atoms(
        Vector{Int64}(topology.angleIndList[:,1]),
        Vector{Int64}(topology.angleIndList[:,2]),
        Vector{Int64}(topology.angleIndList[:,3]),
        Interaction

    )
    # angles = InteractionList3Atoms(
    #     Vector{Int64}(topology.angleIndList[:,1]),
    #     Vector{Int64}(topology.angleIndList[:,2]),
    #     Vector{Int64}(topology.angleIndList[:,3]),
    #     [HarmonicAngle(k=700.0u"kJ * mol^-1", θ0=topology.angleIndList[j,4]*pi/180.0) for j in 1:n_angles_tot]

    # )

    return angles
end


function makeDihedrals(atoms,topology,atomTypeArray,TorsionTypes)
    n_atoms = size(atoms,1)
    n_dihedral_tot = size(topology.dihedralIndList,1)
    
    torsionLog = []
    

    Interaction = []
    for j in 1:n_dihedral_tot
        id1 = Int64(topology.dihedralIndList[j,1])
        id2 = Int64(topology.dihedralIndList[j,2])
        id3 = Int64(topology.dihedralIndList[j,3])
        id4 = Int64(topology.dihedralIndList[j,4])

        element1 = atomTypeArray[id1]["element"]
        hyb1 = atomTypeArray[id1]["hybridization"]
        class1 = atomTypeArray[id1]["class"]
        gafftype1 = atomTypeArray[id1]["gafftype"]

        element2 = atomTypeArray[id2]["element"]
        hyb2 = atomTypeArray[id2]["hybridization"]
        class2 = atomTypeArray[id2]["class"]
        gafftype2 = atomTypeArray[id2]["gafftype"]

        element3 = atomTypeArray[id3]["element"]
        hyb3 = atomTypeArray[id3]["hybridization"]
        class3 = atomTypeArray[id3]["class"]
        gafftype3 = atomTypeArray[id3]["gafftype"]

        element4 = atomTypeArray[id4]["element"]
        hyb4 = atomTypeArray[id4]["hybridization"]
        class4 = atomTypeArray[id4]["class"]
        gafftype4 = atomTypeArray[id4]["gafftype"]

        #For dihedrals I want to find type1, type2, type3, type4 OR type4, type3, type2, type1
        selTorsionTypes = filter(d -> (d["type1"] == gafftype1 && d["type2"] == gafftype2 && d["type3"] == gafftype3 && d["type4"] == gafftype4) ||
                                      (d["type1"] == gafftype4 && d["type2"] == gafftype3 && d["type3"] == gafftype2 && d["type4"] == gafftype1), 
                                      TorsionTypes)
        #but I also need to contemplate wildcards, as not all four atoms may be defined. In that case
        #I want to log some error, so I know which sets of four atoms are not defined
        
        if isempty(selTorsionTypes)
            selTorsionTypes = filter(d -> (d["type1"] == "" && d["type2"] == gafftype2 && d["type3"] == gafftype3 && d["type4"] == "") ||
                                          (d["type1"] == "" && d["type2"] == gafftype3 && d["type3"] == gafftype2 && d["type4"] == ""), 
                                          TorsionTypes)
            msg = "Dihedral parameters for atoms $id1 ($gafftype1), $id2 ($gafftype2), $id3 ($gafftype3) and $id4 ($gafftype4) were not found in AngleTypes. Wildcards were used."
            push!(torsionLog,msg)
        end

        if !isempty(selTorsionTypes)
            #27sep2026: I need an if-statement that checks if more than one periodicity exists
            #Something like if(contains(selTorsionTypes,"periodicity2")
            #do this, etc

            #This is not fixed yet
            if haskey(selTorsionTypes[1],"periodicity1")
                if haskey(selTorsionTypes[1],"periodicity2")
                    push!(Interaction,PeriodicTorsion(periodicities = [Int(selTorsionTypes[1]["periodicity1"]),Int(selTorsionTypes[1]["periodicity2"])], phases = [(selTorsionTypes[1]["phase1"])u"rad",(selTorsionTypes[1]["phase2"])u"rad"], ks=[(selTorsionTypes[1]["k1"])u"kJ/mol",(selTorsionTypes[1]["k2"])u"kJ/mol"]))
                else
                    push!(Interaction,PeriodicTorsion(periodicities = [Int(selTorsionTypes[1]["periodicity1"])], phases = [(selTorsionTypes[1]["phase1"])u"rad"], ks=[(selTorsionTypes[1]["k1"])u"kJ/mol"]))
                end
            end
            
        else
            errorMsg = "Dihedral parameters for atoms $gafftype1, $gafftype2, $gafftype3 and $gafftype4 were not found in AngleTypes."
            push!(torsionLog,errorMsg)
        end
    end

    torsions = InteractionList4Atoms(
        Vector{Int64}(topology.dihedralIndList[:,1]),
        Vector{Int64}(topology.dihedralIndList[:,2]),
        Vector{Int64}(topology.dihedralIndList[:,3]),
        Vector{Int64}(topology.dihedralIndList[:,4]),
        Interaction

    )


    # torsions = InteractionList4Atoms(
    #     Vector{Int64}(topology.dihedralIndList[:,1]),
    #     Vector{Int64}(topology.dihedralIndList[:,2]),
    #     Vector{Int64}(topology.dihedralIndList[:,3]),
    #     Vector{Int64}(topology.dihedralIndList[:,4]),
    #     #PeriodicTorsion(; periodicities, phases, ks, proper)
    #     #I test if the previous separate definition of interaction works instead of the following code:
    #     #[PeriodicTorsion(periodicities = [1], phases = [topology.dihedralIndList[j,5]*pi/180.0], ks=[7000.0u"kJ * mol^-1"]) for j in 1:n_dihedral_tot]
    #     Interaction
    # )

    return torsions, torsionLog
end

function determineAtomTypes(atoms,topology)
    n_atoms = size(atoms,1)
    #need to decide if I want a dictionary or an Array with Any types
    #If a Dict only accepts a key and value, I may need to nest Dicts
    #atomTypeDict = Dict()
    #I think the solution is to use an array of Dicts
    #Each Dict must tell which element it is, and which hybridization and class it belongs to (sp2, sp3, carbonyl, -CH2-, -C=O-, etc)
    atomTypeArray = Any[]
    #my_Array = Any[]
    #push!(my_array, [1, 2, 3])       # Adds a Vector{Int64}
    #push!(my_array, ["a", "b"])      # Adds a Vector{String}
    #push!(my_array, 42)              # Adds an Int64

    for j in 1:n_atoms
        mass = ustrip(atoms[j].mass) #this assumes it is given in g mol^-1

        element = keys(filter(p -> p.second[2] == mass, elementTypeDict))
        if first(element) == "H"
            atomTypeDict = Dict(
                "element" => "H",
                "hybridization" => "s",
                "class" => "H",
                "gafftype" => "",
                )
            number_neighbors, indices, neighbor_elements, neighbor_indices = findNeighbors(j,atoms,topology,elementTypeDict)

            if number_neighbors["C"] == 1
                atomTypeDict["class"] = "H-C"
                atomTypeDict["gafftype"] = "hc" #on aliphatic carbon, so far no distinction with aromatic carbon
            elseif number_neighbors["N"] == 1
                atomTypeDict["class"] = "H-N"
                atomTypeDict["gafftype"] = "hn"
            elseif number_neighbors["O"] == 1
                atomTypeDict["class"] = "H-O"
                atomTypeDict["gafftype"] = "ho"
            end

            push!(atomTypeArray,atomTypeDict)
        elseif first(element) == "C"
            atomTypeDict = Dict(
                "element" => "C",
                "hybridization" => "",
                "class" => "",
                "gafftype" => "",
                )
            
            number_neighbors, indices, neighbor_elements, neighbor_indices = findNeighbors(j,atoms,topology,elementTypeDict)
            # number_H_neighbors = number_neighbors["H"]
            # number_C_neighbors = number_neighbors["C"] 
            # number_O_neighbors = number_neighbors["O"] 
            # println(number_neighbors["H"])
            # println(number_neighbors["C"])
            # println(number_neighbors["O"])


            if length(indices) == 3
                angleRowInd,angleValues = getAngles(j,neighbor_indices,topology)
                # println("Output from getAngles for:")
                # println(j)
                # println(neighbor_indices)
                # println(angleRowInd)
                # println(angleValues)
                if sum(angleValues) > 350.0 #for a planar site they should add to exactly 360 degrees
                    atomTypeDict["hybridization"] = "sp2"
                else
                    atomTypeDict["hybridization"] = "sp3*" #the star is a placeholder, it could mean an empty site from a missing atom
                end

                if number_neighbors["H"] == 2
                    atomTypeDict["class"] = "=CH2"
                    atomTypeDict["gafftype"] = "c2"
                elseif number_neighbors["O"] == 2
                    atomTypeDict["class"] = "-CO2-"
                    atomTypeDict["gafftype"] = "c" #gaff carbonyl carbon
                elseif number_neighbors["C"] == 2
                    atomTypeDict["class"] = "-CX="
                    atomTypeDict["gafftype"] = "ca" #gaff aromatic carbon (just a guess so far)
                elseif number_neighbors["C"] == 3
                    atomTypeDict["class"] = "-C(bridge)="
                    atomTypeDict["gafftype"] = "cp" #gaff aromatic carbon, bridging phenyl units (just a guess so far)
                end
            elseif length(indices) == 4
                atomTypeDict["hybridization"] = "sp3"
                if number_neighbors["H"] == 3
                    atomTypeDict["class"] = "CH3-"
                    atomTypeDict["gafftype"] = "c3"
                elseif number_neighbors["O"] == 1
                    atomTypeDict["class"] = "CX2O-"
                    atomTypeDict["gafftype"] = "c3"
                end
            end
            
            push!(atomTypeArray,atomTypeDict)
        elseif first(element) == "N" #to do
            atomTypeDict = Dict(
                "element" => "N",
                "hybridization" => "",
                "class" => "",
                "gafftype" => "",
                )


            push!(atomTypeArray,atomTypeDict)    
        elseif first(element) == "O" #to do
            atomTypeDict = Dict(
                "element" => "O",
                "hybridization" => "",
                "class" => "",
                "gafftype" => "",
                )

            number_neighbors, indices = findNeighbors(j,atoms,topology,elementTypeDict)
            if length(indices) == 1
                atomTypeDict["hybridization"] = "sp2"
                if number_neighbors["C"] == 1
                    atomTypeDict["class"] = "C=O"
                    atomTypeDict["gafftype"] = "o" #gaff carbonyl oxygen
                end
            elseif length(indices) == 2
                atomTypeDict["hybridization"] = "sp3"
                if number_neighbors["C"] == 1 && number_neighbors["H"] == 1
                    atomTypeDict["class"] = "C-O-H"
                    atomTypeDict["gafftype"] = "oh" #gaff alcohol oxygen
                elseif number_neighbors["C"] == 2
                    atomTypeDict["class"] = "C-O-C"
                    atomTypeDict["gafftype"] = "os" #gaff ether or ester oxygen    
                end
            end
            
            push!(atomTypeArray,atomTypeDict)    

        end

    end

    return atomTypeArray


end


function findNeighbors(j,atoms,topology,elementTypeDict)
    indices = findall(x -> x == j, topology.bondIndList)    #j will be cast as Float by itself
    #as the second coordinate in each element of the indices vector can be either 1 or 2, we can find the index for the other atom
    col_neighbors = []
    row_neighbors = []
    for k in 1:length(indices)
        indices[k][2] == 1 ? push!(row_neighbors, 2) : push!(row_neighbors, 1)
        push!(col_neighbors,indices[k][1])
    end

    neighbor_elements = []
    neighbor_indices = []
    for k in 1:length(col_neighbors)
        idx_neighbor = Int64(topology.bondIndList[col_neighbors[k],row_neighbors[k]])
        mass_neighbor = ustrip(atoms[idx_neighbor].mass) #this assumes it is given in g mol^-1
        neighbor_element = keys(filter(p -> p.second[2] == mass_neighbor, elementTypeDict))
        push!(neighbor_elements,first(neighbor_element))
        push!(neighbor_indices,idx_neighbor)
    end
    # println("The neighbor elements for this C are: ")
    # println(neighbor_elements)
    number_neighbors = Dict(
        "H" => 0,
        "C" => 0,
        "N" => 0,
        "O" => 0,
    )
    number_neighbors["H"] = length(findall(x -> x == "H", neighbor_elements))
    number_neighbors["C"] = length(findall(x -> x == "C", neighbor_elements))
    number_neighbors["O"] = length(findall(x -> x == "O", neighbor_elements))
    number_neighbors["N"] = length(findall(x -> x == "N", neighbor_elements))
    
    #I think indices and neighbor_indices are the same
    #actually not the same, indices are the coordinates to find the neighbor_indices in bondIndList
    return number_neighbors, indices, neighbor_elements, neighbor_indices
end

function getAngles(current_index,neighbor_indices,topology)

    #The following line was suggested by Gemini. It seems to work
    #It now assumes that neighbor_indices are given in increasing order, AND angles in angleIndList always have the lower index as the first column
    if length(neighbor_indices) == 2
        angleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[1] current_index neighbor_indices[2]]), 1:size(topology.angleIndList, 1))
        if isempty(angleRowInd)
            angleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[2] current_index neighbor_indices[1]]), 1:size(topology.angleIndList, 1))
        end

        angleValues = topology.angleIndList[angleRowInd,4]
    elseif length(neighbor_indices) == 3
        angleRowInd = []
        angleValues = []

        #1 and #2
        currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[1] current_index neighbor_indices[2]]), 1:size(topology.angleIndList, 1))
        if isempty(angleRowInd)
            currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[2] current_index neighbor_indices[1]]), 1:size(topology.angleIndList, 1))
        end        
        if !isempty(currentAngleRowInd)
            push!(angleRowInd,first(currentAngleRowInd))
            push!(angleValues,first(topology.angleIndList[currentAngleRowInd,4]))
        end
        #1 and #3
        currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[1] current_index neighbor_indices[3]]), 1:size(topology.angleIndList, 1))
        if isempty(angleRowInd)
            currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[3] current_index neighbor_indices[1]]), 1:size(topology.angleIndList, 1))
        end
        if !isempty(currentAngleRowInd)
            push!(angleRowInd,first(currentAngleRowInd))
            push!(angleValues,first(topology.angleIndList[currentAngleRowInd,4]))
        end
        #2 and #3
        currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[2] current_index neighbor_indices[3]]), 1:size(topology.angleIndList, 1))
        if isempty(angleRowInd)
            currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[3] current_index neighbor_indices[2]]), 1:size(topology.angleIndList, 1))
        end
        if !isempty(currentAngleRowInd)
            push!(angleRowInd,first(currentAngleRowInd))
            push!(angleValues,first(topology.angleIndList[currentAngleRowInd,4]))
        end

    elseif length(neighbor_indices) == 4
        angleRowInd = []
        angleValues = []

        #1 and #2
        currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[1] current_index neighbor_indices[2]]), 1:size(topology.angleIndList, 1))
        if isempty(angleRowInd)
            currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[2] current_index neighbor_indices[1]]), 1:size(topology.angleIndList, 1))
        end
        if !isempty(currentAngleRowInd)
            push!(angleRowInd,first(currentAngleRowInd))
            push!(angleValues,first(topology.angleIndList[currentAngleRowInd,4]))
        end
        #1 and #3
        currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[1] current_index neighbor_indices[3]]), 1:size(topology.angleIndList, 1))
        if isempty(angleRowInd)
            currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[3] current_index neighbor_indices[1]]), 1:size(topology.angleIndList, 1))
        end
        if !isempty(currentAngleRowInd)
            push!(angleRowInd,first(currentAngleRowInd))
            push!(angleValues,first(topology.angleIndList[currentAngleRowInd,4]))
        end
        #1 and #4
        currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[1] current_index neighbor_indices[4]]), 1:size(topology.angleIndList, 1))
        if isempty(angleRowInd)
            currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[4] current_index neighbor_indices[1]]), 1:size(topology.angleIndList, 1))
        end
        if !isempty(currentAngleRowInd)
            push!(angleRowInd,first(currentAngleRowInd))
            push!(angleValues,first(topology.angleIndList[currentAngleRowInd,4]))
        end
        #2 and #3
        currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[2] current_index neighbor_indices[3]]), 1:size(topology.angleIndList, 1))
        if isempty(angleRowInd)
            currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[2] current_index neighbor_indices[2]]), 1:size(topology.angleIndList, 1))
        end
        if !isempty(currentAngleRowInd)
            push!(angleRowInd,first(currentAngleRowInd))
            push!(angleValues,first(topology.angleIndList[currentAngleRowInd,4]))
        end  
        #2 and #4
        currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[2] current_index neighbor_indices[4]]), 1:size(topology.angleIndList, 1))
        if isempty(angleRowInd)
            currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[4] current_index neighbor_indices[2]]), 1:size(topology.angleIndList, 1))
        end
        if !isempty(currentAngleRowInd)
            push!(angleRowInd,first(currentAngleRowInd))
            push!(angleValues,first(topology.angleIndList[currentAngleRowInd,4]))
        end
        #3 and #4
        currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[3] current_index neighbor_indices[4]]), 1:size(topology.angleIndList, 1))
        if isempty(angleRowInd)
            currentAngleRowInd = findall(i -> @view(topology.angleIndList[i, 1:3]) == vec([neighbor_indices[4] current_index neighbor_indices[3]]), 1:size(topology.angleIndList, 1))
        end
        if !isempty(currentAngleRowInd)
            push!(angleRowInd,first(currentAngleRowInd))
            push!(angleValues,first(topology.angleIndList[currentAngleRowInd,4]))
        end                    
    end

    return angleRowInd,angleValues
end

#Energy calculations################################################################################
function calcCoulomb(frames,atomPositions,atomNames,atomCharges)
    

end


#Output and File Formats############################################################################
function printXYZ(fragment::SimpleFragment)

    nAtoms = fragment.nAtoms

    for k = 1:nAtoms
        mass = fragment.Masses[k]
        element = keys(filter(p -> p.second[2] == mass, elementTypeDict))
        xk = round(fragment.Coords[k,1]+0e-6, digits = 6)
        yk = round(fragment.Coords[k,2]+0e-6, digits = 6)
        zk = round(fragment.Coords[k,3]+0e-6, digits = 6)

        @printf("%s   %f   %f   %f \n", first(element), xk, yk, zk)
        #println(first(element), "\t", xk, "\t", yk, "\t", zk)
    end


end

function printPDB(fragment::SimpleFragment, )


end

#This function is meant to read and parse forcefield xml files using EzXML
function readFFFile(filename::String)
    ff = readxml(filename)
    ForceField = root(ff)

    fields = []
    AtomTypes = []
    BondTypes = []
    AngleTypes = []
    TorsionTypes = []
    NonBondedTypes = []
    for field in eachelement(ForceField)
        #field_name = field["name"]
        #field_name = nodecontent(field)
        push!(fields,field)
        
            if nodename(field) == "AtomTypes"
                
                for atomtype in eachelement(field)
                    #println(atomtype)
                    push!(AtomTypes,atomtype)
                end
                #println(field)
            elseif nodename(field) == "HarmonicBondForce"
                for bondtype in eachelement(field)
                    #example: <Bond type1="n2" type2="sy" length="0.1544" k="411220.256"/>
                    bondtype_content = split(string(bondtype))
                    #println(bondtype_content)
                    if length(bondtype_content) == 5
                        type1_field = bondtype_content[2]
                        type1_field = replace(type1_field, r"[\"\\/]" => "")
                        type1 = split(type1_field,"=")
                        #println(type1[2])

                        type2_field = bondtype_content[3]
                        type2_field = replace(type2_field, r"[\"\\/]" => "")
                        type2 = split(type2_field,"=")

                        length_field = bondtype_content[4]
                        length_field = replace(length_field, r"[\"\\/]" => "")
                        bondlength = split(length_field,"=")
                        #println(bondlength)
                        bondlength_value = parse(Float64,bondlength[2])

                        k_field = bondtype_content[5]
                        k_field = replace(k_field, r"[\"\\/>]" => "")
                        bond_k = split(k_field,"=")
                        #println(bond_k)
                        k_value = parse(Float64,bond_k[2])
                    end
                    bondType = Dict(
                        "type1" => type1[2],
                        "type2" => type2[2],
                        "length" => bondlength_value,
                        "k" => k_value
                    )

                    #println(bondtype_content)
                    push!(BondTypes,bondType)
                end
            elseif nodename(field) == "HarmonicAngleForce"
                for angletype in eachelement(field)
                    #example: <Angle type1="c2" type2="c2" type3="s4" angle="2.09160257559" k="541.802896"/>
                    
                    angletype_content = split(string(angletype))
                    #println(bondtype_content)
                    if length(angletype_content) == 6
                        type1_field = angletype_content[2]
                        type1_field = replace(type1_field, r"[\"\\/]" => "")
                        type1 = split(type1_field,"=")
                        #println(type1[2])

                        type2_field = angletype_content[3]
                        type2_field = replace(type2_field, r"[\"\\/]" => "")
                        type2 = split(type2_field,"=")

                        type3_field = angletype_content[4]
                        type3_field = replace(type3_field, r"[\"\\/]" => "")
                        type3 = split(type3_field,"=")

                        angle_field = angletype_content[5]
                        angle_field = replace(angle_field, r"[\"\\/]" => "")
                        bondangle = split(angle_field,"=")
                        #println(bondlength)
                        angle_value = parse(Float64,bondangle[2])

                        k_field = angletype_content[6]
                        k_field = replace(k_field, r"[\"\\/>]" => "")
                        angle_k = split(k_field,"=")
                        #println(bond_k)
                        k_value = parse(Float64,angle_k[2])
                    end
                    angleType = Dict(
                        "type1" => type1[2],
                        "type2" => type2[2],
                        "type3" => type3[2],
                        "angle" => angle_value,
                        "k" => k_value
                    )

                    #println(bondtype_content)
                    push!(AngleTypes,angleType)
                end                
            elseif nodename(field) == "PeriodicTorsionForce"
                for torsiontype in eachelement(field)
                    #example: <Proper type1="c3" type2="p3" type3="c" type4="c3" periodicity1="2" phase1="3.14159265359" k1="6.434992"/>
                    #example: <Improper type1="ca" type2="" type3="" type4="ha" periodicity1="2" phase1="3.14159265359" k1="4.6024"/>
                    
                    torsiontype_content = split(string(torsiontype))
                    #println(bondtype_content)
                    if length(torsiontype_content) >= 8 #but there can be more
                        type1_field = torsiontype_content[2]
                        type1_field = replace(type1_field, r"[\"\\/]" => "")
                        type1 = split(type1_field,"=")
                        #println(type1[2])

                        type2_field = torsiontype_content[3]
                        type2_field = replace(type2_field, r"[\"\\/]" => "")
                        type2 = split(type2_field,"=")

                        type3_field = torsiontype_content[4]
                        type3_field = replace(type3_field, r"[\"\\/]" => "")
                        type3 = split(type3_field,"=")

                        type4_field = torsiontype_content[5]
                        type4_field = replace(type4_field, r"[\"\\/]" => "")
                        type4 = split(type4_field,"=")

                        period1_field = torsiontype_content[6]
                        period1_field = replace(period1_field, r"[\"\\/]" => "")
                        period1 = split(period1_field,"=")
                        
                        period1_value = parse(Float64,period1[2])

                        phase1_field = torsiontype_content[7]
                        phase1_field = replace(phase1_field, r"[\"\\/]" => "")
                        phase1 = split(phase1_field,"=")
                        
                        phase1_value = parse(Float64,phase1[2])


                        k1_field = torsiontype_content[8]
                        k1_field = replace(k1_field, r"[\"\\/>]" => "")
                        torsion_k1 = split(k1_field,"=")
                        #println(bond_k)
                        k1_value = parse(Float64,torsion_k1[2])

                        if length(torsiontype_content) >= 11
                            period2_field = torsiontype_content[9]
                            period2_field = replace(period2_field, r"[\"\\/]" => "")
                            period2 = split(period2_field,"=")                     
                            period2_value = parse(Float64,period2[2])

                            phase2_field = torsiontype_content[10]
                            phase2_field = replace(phase2_field, r"[\"\\/]" => "")
                            phase2 = split(phase2_field,"=")
                            phase2_value = parse(Float64,phase2[2])

                            k2_field = torsiontype_content[11]
                            k2_field = replace(k2_field, r"[\"\\/>]" => "")
                            torsion_k2 = split(k2_field,"=")
                            k2_value = parse(Float64,torsion_k2[2])
                        end
                    end
                    if @isdefined period2_value
                        torsionType = Dict(
                            "type1" => type1[2],
                            "type2" => type2[2],
                            "type3" => type3[2],
                            "type4" => type4[2],
                            "periodicity1" => period1_value,
                            "phase1" => phase1_value,
                            "k1" => k1_value,
                            "periodicity2" => period2_value,
                            "phase2" => phase2_value,
                            "k2" => k2_value
                        )
                    else
                        torsionType = Dict(
                            "type1" => type1[2],
                            "type2" => type2[2],
                            "type3" => type3[2],
                            "type4" => type4[2],
                            "periodicity1" => period1_value,
                            "phase1" => phase1_value,
                            "k1" => k1_value
                        )
                    end
                    #println(bondtype_content)
                    push!(TorsionTypes,torsionType)
                end
            elseif nodename(field) == "NonbondedForce"
                for nonbondedtype in eachelement(field)
                    nonbondedtype_content = split(string(nonbondedtype))
                    #println(nonbondedtype_content)
                    if length(nonbondedtype_content) == 4 
                        typeLJ_field = nonbondedtype_content[2]
                        typeLJ_field = replace(typeLJ_field, r"[\"\\/]" => "")
                        typeLJ = split(typeLJ_field,"=")

                        sigma_field = nonbondedtype_content[3]
                        sigma_field = replace(sigma_field, r"[\"\\/]" => "")
                        sigma = split(sigma_field,"=")
                        sigma_value = parse(Float64,sigma[2])

                        epsilon_field = nonbondedtype_content[4]
                        epsilon_field = replace(epsilon_field, r"[\"\\/>]" => "")
                        epsilon = split(epsilon_field,"=")
                        epsilon_value = parse(Float64,epsilon[2])
                    
                        nonbondedType = Dict(
                            "typeLJ" => typeLJ[2],
                            "sigma" => sigma_value,
                            "epsilon" => epsilon_value
                        )
                        push!(NonBondedTypes,nonbondedType)
                    end #the definition of the Dict is inside the if-statement. This is different than for the other fields.
                end
            end

        #end
        #println(field)
    end
    return ff, fields, AtomTypes, BondTypes, AngleTypes, TorsionTypes, NonBondedTypes
end

#ff = MolecularForceField(Float32,"gaff2.xml")