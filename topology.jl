include("AtomParameters.jl")

mutable struct MoleculeTopology
    nAtoms::Int
    atomIndList::Matrix{Any}
    bondIndList::Union{Matrix{Int64},Matrix{Float64}}
    angleIndList::Matrix{Any}
    dihedralIndList::Matrix{Any}
    improperIndList
    nonBondedIndList::Matrix{Any}

end

function buildTopology(molecule,elementTypeDict,atomRadiiDict)
    #From now on, distances will be in nm
    #atomRadiiDict are given in pm, so they need to be divided by 1000 to obtain nm
    nAtoms = size(molecule.Masses,1)
    #This will be removed
    bondLengthThreshold = 0.16 #This should be later changed to account for atom types

    
    atomIndList = zeros(nAtoms,3)
    bondIndList = zeros(nAtoms*(nAtoms-1),4)
    angleIndList = zeros(nAtoms*(nAtoms-1)*(nAtoms-2),5)
    dihedralIndList = zeros(nAtoms*(nAtoms-1)*(nAtoms-2)*(nAtoms-3),6)
    
    
    for k in 1:nAtoms
        #atomIndList = [atomicNumber sigma epsilon]

        moleculeNumbers = floor.(Int64,molecule.Masses/2 .+ 0.1) #This is a hack for finding the atom numbers from the masses. It will not work when many heavy elements are included
        moleculeNumbers[moleculeNumbers.==0].=1 #This fixes H atoms, which otherwise would have zero atomic number
        
        #This numbers have to be changed to element specific numbers
        σ_List = [1. 1. 1. 1. 1. 1. 1. 1. 1. 1. ]
        ϵ_List = [1. 1. 1. 1. 1. 1. 1. 1. 1. 1. ] #Which units? Have to be consistent overall. Better to use kJ/mol

        atomNumber = Int64(moleculeNumbers[k])

        atomIndList[k,:] = [atomNumber σ_List[atomNumber] ϵ_List[atomNumber]]
    end

    currentBondMarker = 1
    for k in 1:nAtoms
        atom1 = molecule.Coords[k,:]
        #This works, after trying a few things (the two is the second element in the vector [atomNumber, atomMass, sigma, epsilon])
        element1 = keys(filter(p -> p.second[2] == molecule.Masses[k], elementTypeDict))
        atomRadius1 = atomRadiiDict[first(element1)]


    bondLengthThreshold = 0.16 #This should be later changed to account for atom types

        for l in k+1:nAtoms
            atom2 = molecule.Coords[l,:]
            element2 = keys(filter(p -> p.second[2] == molecule.Masses[l], elementTypeDict))
            atomRadius2 = atomRadiiDict[first(element2)]

            bondLengthThreshold1_2 = (first(atomRadius1) + first(atomRadius2))/1000 #pm to nm
            bond1_2 = atom2-atom1
            # println("current bond length")
            # println(norm(bond1_2))
            # println("threshold for pair of elements")
            # println(bondLengthThreshold1_2)
            
            #These are dummy variables which later have to be changed by forcefield parameters
            minDistance = norm(bond1_2)
            bondForceConstant = 0.0001

            if norm(bond1_2) < bondLengthThreshold1_2
                bondIndList[currentBondMarker,:] = [k l minDistance bondForceConstant]
                currentBondMarker += 1
            end
        end
    end

    bondIndList = bondIndList[1:currentBondMarker-1,:]
    nBonds = size(bondIndList,1)
    
    currentAngleMarker = 1
    for kl in 1:nBonds
        id1 = Int64(bondIndList[kl,1])
        id2 = Int64(bondIndList[kl,2])

        atom1 = molecule.Coords[id1,:]
        atom2 = molecule.Coords[id2,:]  #id2 is always supposed to be larger than id1

        element1 = keys(filter(p -> p.second[2] == molecule.Masses[id1], elementTypeDict))
        atomRadius1 = atomRadiiDict[first(element1)]

        element2 = keys(filter(p -> p.second[2] == molecule.Masses[id2], elementTypeDict))
        atomRadius2 = atomRadiiDict[first(element2)]

        bondLengthThreshold1_2 = (first(atomRadius1) + first(atomRadius2))/1000 #pm to nm

        #but what if atom1 is bonded to a third atom??

        bond2_1 = atom1-atom2
        n_bond2_1 = bond2_1/norm(bond2_1) #normalized vector in the 2->1 direction


        for m in 1:nAtoms #I think this also needs to run from 1 to nAtoms, but check that m neq id1 or id2
            
            if (m == id1) || (m == id2)
                continue
            end

            atom3 = molecule.Coords[m,:]
            element3 = keys(filter(p -> p.second[2] == molecule.Masses[m], elementTypeDict))
            atomRadius3 = atomRadiiDict[first(element3)]

            bondLengthThreshold1_3 = (first(atomRadius1) + first(atomRadius3))/1000 #pm to nm
            bondLengthThreshold2_3 = (first(atomRadius2) + first(atomRadius3))/1000 #pm to nm

            bond2_3 = atom3-atom2
            n_bond2_3 = bond2_3/norm(bond2_3) #normalized vector in the 2->3 direction
            currentAngle1_2_3 = acos(dot(n_bond2_1,n_bond2_3))*180/pi #angle from dot product of two vectors (2->1 dot d->3)

            bond1_3 = atom3-atom1
            n_bond1_3 = bond1_3/norm(bond1_3) #normalized vector in the 2->3 direction
            #here I reverse the sign of the n_bond2_1, as now atom1 is the center
            currentAngle2_1_3 = acos(dot(-n_bond2_1,n_bond1_3))*180/pi #angle from dot product of two vectors (2->1 dot d->3)

            angleForceConstant = 0.00001
            
            if norm(bond2_3) < bondLengthThreshold2_3
                angleIndList[currentAngleMarker,:] = [id1 id2 m currentAngle1_2_3 angleForceConstant]
                currentAngleMarker += 1
            elseif norm(bond1_3) < bondLengthThreshold1_3
                angleIndList[currentAngleMarker,:] = [id2 id1 m currentAngle2_1_3 angleForceConstant]
                currentAngleMarker += 1
            end
        end

    end

    angleIndList = angleIndList[1:currentAngleMarker-1,:]

    currentDihedralMarker = 1
    for kl in 1:nBonds
        id1 = Int64(bondIndList[kl,1])
        id2 = Int64(bondIndList[kl,2])

        atom1 = molecule.Coords[id1,:]
        atom2 = molecule.Coords[id2,:]  #id2 is always supposed to be larger than id1

        bond2_1 = atom1-atom2
        n_bond2_1 = bond2_1/norm(bond2_1)

        for mn in kl+1:nBonds
            id3 = Int64(bondIndList[mn,1])
            id4 = Int64(bondIndList[mn,2])
            
            atom3 = molecule.Coords[id3,:]
            atom4 = molecule.Coords[id4,:]

            bond2_3 = atom3-atom2
            
            #bond2_3 = atom3-atom2
            n_bond2_3 = bond2_3/norm(bond2_3)
            
            bond2_4 = atom4-atom2 #atom2 may be bonded to atom4, not atom3, we do not know, depends on the numbering of atoms, which for a branched molecule may involve numbers which don't always increment from left to right
            n_bond2_4 = bond2_4/norm(bond2_4)
            #currentAngle = acos(dot(n_bond2_1,n_bond2_3))*180/pi
            
            dihedralForceConstant = 0.00001
            if (id2 != id3) && (id2 != id4) && (id1 != id3) && (id1 != id4)
                if (norm(bond2_3) < bondLengthThreshold)  #if there is a bond between 2 and 3 (but they are not the same atom)
                    normal123 = cross(n_bond2_1,n_bond2_3)
                    n_normal_123 = normal123/norm(normal123)

                    bond3_4 = atom4-atom3
                    n_bond3_4 = bond3_4/norm(bond3_4)

                    normal234 = cross(-n_bond2_3,n_bond3_4)
                    n_normal_234 = normal234/norm(normal234)

                    # println("Debugging dihedral:")
                    # println([id1 id2 id3 id4])
                    # println([n_normal_123 n_normal_234])
                    
                    dot_123_234 = dot(n_normal_123,n_normal_234)
                    if abs(dot_123_234) > 1.0
                        currentDihedral = 180.0
                    else
                        currentDihedral = acos(dot_123_234)*180/pi
                    end

                    dihedralIndList[currentDihedralMarker,:] = [id1 id2 id3 id4 currentDihedral dihedralForceConstant]
                    currentDihedralMarker += 1

                elseif (norm(bond2_4) < bondLengthThreshold)  #if there is a bond between 2 and 3 (but they are not the same atom)
                    normal124 = cross(n_bond2_1,n_bond2_4)
                    n_normal_124 = normal124/norm(normal124)

                    bond4_3 = atom3-atom4
                    n_bond4_3 = bond4_3/norm(bond4_3)

                    normal243 = cross(-n_bond2_3,n_bond4_3)
                    n_normal_243 = normal243/norm(normal243)

                    # println("Debugging dihedral:")
                    # println([id1 id2 id3 id4])
                    # println([n_normal_123 n_normal_234])
                    
                    dot_124_243 = dot(n_normal_124,n_normal_243)
                    if abs(dot_124_243) > 1.0
                        currentDihedral = 180.0
                    else
                        currentDihedral = acos(dot_124_243)*180/pi
                    end

                    dihedralIndList[currentDihedralMarker,:] = [id1 id2 id4 id3 currentDihedral dihedralForceConstant]
                    currentDihedralMarker += 1
                end
            end
        end

    end

    dihedralIndList = dihedralIndList[1:currentDihedralMarker-1,:]
    
    improperIndList = []
    #This code is a simplified version of the code in findNeighbors
    #I cannot call findNeighbors, because one of its arguments is topology (but only bondIndList is needed)
    for k in 1:nAtoms
        indices = findall(x -> x == k, bondIndList)
        col_neighbors = []
        row_neighbors = []
        for k in 1:length(indices)
            indices[k][2] == 1 ? push!(row_neighbors, 2) : push!(row_neighbors, 1)
            push!(col_neighbors,indices[k][1])
        end
        neighbor_indices = []

        for k in 1:length(col_neighbors)
            idx_neighbor = Int64(bondIndList[col_neighbors[k],row_neighbors[k]])
            push!(neighbor_indices,idx_neighbor)
        end

        if length(neighbor_indices) == 3 #a central atom k, bound to 3 atoms
            push!(improperIndList,[neighbor_indices[1] k neighbor_indices[2] neighbor_indices[3]])
        end

    end
    
    topology = MoleculeTopology(
        nAtoms,
        atomIndList,
        bondIndList,
        angleIndList,
        dihedralIndList,
        improperIndList,
        zeros(1,4),
    )
    

    return topology
end

function buildFragmentConnectivity(topology::MoleculeTopology,atom1::Int64,atom2::Int64)
#This function separates a molecule in two parts, by breaking a bond connecting atom1 and atom2, and returs two lists of indices of the atoms of each fragment
#the lists should be exclusive (no repeated atoms, and all atoms accounted for)
#It has to search all bonds to determine the two fragments of a molecule. If the broken bond does not separate the molecule in two parts,
#it should tell in some way.
#Example:
# topology.bondIndList
# 17×4 Matrix{Any}:
#   1.0   3.0  1.44     0.0001
#   2.0   9.0  1.21108  0.0001
#   3.0   4.0  1.09698  0.0001
#   3.0   5.0  1.52631  0.0001
#   3.0   9.0  1.51306  0.0001
#   5.0   6.0  1.09356  0.0001
#   5.0   7.0  1.09296  0.0001
#   5.0   8.0  1.09296  0.0001
#   9.0  10.0  1.35756  0.0001
#  10.0  12.0  1.44     0.0001
#  11.0  18.0  1.21108  0.0001
#  12.0  13.0  1.09698  0.0001
#  12.0  14.0  1.52631  0.0001
#  12.0  18.0  1.51306  0.0001
#  14.0  15.0  1.09356  0.0001
#  14.0  16.0  1.09296  0.0001
#  14.0  17.0  1.09296  0.0001

#I want the indices ordered, so that atom2 is larger
if atom1 > atom2
    atom1, atom2 = atom2, atom1
end

bondIndList = topology.bondIndList[:,1:2]
next1 = Float64(atom1)
next2 = Float64(atom2)

#This whole function is unnecessarily complicated.
#I can move the rows in bondIndList, that contain atom2 in the second position,
#right after the row which contains [atom1 atom2] (Let's say row splitInd)
#Then I just split the bondIndList into bondIndList[1:splitInd-1,:] and bondIndList[splitInd+1:end,:]
#Then I collapse the 2d lists into 1-D, and get only the unique indices.
#Then I check that the number of atoms in both lists add up to the original number of atoms

rowInd  = findall(all(bondIndList .== [atom1 atom2], dims=2))

rowNumber = rowInd[1][1] #to get the actual row number for the CartesianIndex which is rowInd

rowsAtom2_is_second = findall(x->x==atom2,bondIndList[1:rowNumber-1,2])

if ~isempty(rowsAtom2_is_second)

    rowContents_Atom2_is_second = bondIndList[rowsAtom2_is_second,:]
    newBondIndList1_temp = bondIndList[1:rowNumber-1,:]
    
    #the way I'm doing this only can work if rowsAtom2_is_second only has one element
    #What I need is a way to remove several rows(indexed by rowsAtom2_is_second) from newBondIndList1_temp
    newBondIndList1 = vcat(newBondIndList1_temp[1:first(rowsAtom2_is_second)-1,:],newBondIndList1_temp[first(rowsAtom2_is_second)+1:end,:])
    
    newBondIndList2 = bondIndList[rowNumber+1:end,:]
    newBondIndList2 = vcat(rowContents_Atom2_is_second,newBondIndList2)
else
    newBondIndList1 = bondIndList[1:rowNumber-1,:]
    newBondIndList2 = bondIndList[rowNumber+1:end,:]
end

indList1 = unique(reshape(newBondIndList1,(length(newBondIndList1),1)))
indList2 = unique(reshape(newBondIndList2,(length(newBondIndList2),1)))

if (length(indList1)+length(indList2)) > topology.nAtoms
    print("The chosen bond breakage does not separate the molecule in two fragments")
end

return bondIndList, indList1, indList2

end

function findBondedAtoms(molecule,elementTypeDict,atomRadiiDict,units::String)
    #This function is abstracted from topology, but is only focused on finding bonds
    nAtoms = size(molecule.Masses,1)

    if units == "nm"
        scaleFactor = 1
    elseif units == "angs"
        scaleFactor = 0.1
    else
        error("Distance units not correctly specified, please use either 'nm' or 'angs'")
    end
    
    bondIndList = zeros(nAtoms*(nAtoms-1),4)
    
    currentBondMarker = 1
    for k in 1:nAtoms
        atom1 = molecule.Coords[k,:].*scaleFactor
        #This works, after trying a few things (the two is the second element in the vector [atomNumber, atomMass, sigma, epsilon])
        element1 = keys(filter(p -> p.second[2] == molecule.Masses[k], elementTypeDict))
        atomRadius1 = atomRadiiDict[first(element1)]


        for l in k+1:nAtoms
            atom2 = molecule.Coords[l,:].*scaleFactor
            element2 = keys(filter(p -> p.second[2] == molecule.Masses[l], elementTypeDict))
            atomRadius2 = atomRadiiDict[first(element2)]

            bondLengthThreshold1_2 = (first(atomRadius1) + first(atomRadius2))/1000 #pm to nm
            bond1_2 = atom2-atom1

            
            #These are dummy variables which later have to be changed by forcefield parameters
            minDistance = norm(bond1_2)
            bondForceConstant = 0.0001

            if norm(bond1_2) < bondLengthThreshold1_2
                bondIndList[currentBondMarker,:] = [k l minDistance bondForceConstant]
                currentBondMarker += 1
            end
        end
    end

    bondIndList = bondIndList[1:currentBondMarker-1,:]
    nBonds = size(bondIndList,1)
    return bondIndList, nBonds
end
