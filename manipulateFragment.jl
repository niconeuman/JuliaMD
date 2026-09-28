using LinearAlgebra



function rodriguesRotations(fragment::Fragment,pivotAtomID::Int64,vectorAxis::Vector{Float64},theta::Float64)

    #This function rotates a matrix of nAtom-by-three numbers around vectorAxis, an angle theta
    #the matrix is obtained from a fragment.Coords field, and used to overwrite this field
    #it uses Rodrigues' rotation formula, given below
    #v_rot = v*cos(theta)+cross(k,v)*sin(theta)+k*(dot(k,v))*(1-cos(theta))
    prevCoordinates = fragment.Coords
    pivotCoords = fragment.Coords[pivotAtomID,:]'
    prevCoordinates = prevCoordinates .- pivotCoords #This is meant to rotate around the origin
    newCoordinates = similar(prevCoordinates)
    nAtoms = size(prevCoordinates,1)

    k = vectorAxis

    #angles are in radians

    for kAtom in 1:nAtoms
        v = prevCoordinates[kAtom,:]
        newCoordinates[kAtom,:] = v*cos(theta)+cross(k,v)*sin(theta)+k*(dot(k,v))*(1-cos(theta))
    end

    fragment.Coords = newCoordinates .+ pivotCoords
end

function eulerRotation(fragment::SimpleFragment,pivotAtomID::Int64,alpha::Float64,beta::Float64,gamma::Float64)

    #This function rotates a matrix of nAtom-by-three numbers using Euler Matrices
    #It first rotates around z, then around x and then again around z.
    #Angles must be given in radian

    # Elementary Rotation Matrices
RotZ_a =  [cos(alpha) -sin(alpha) 0.0;
            sin(alpha)  cos(alpha) 0.0;
            0.0     0.0    1.0] 

RotX_b = [cos(beta)  0.0  sin(beta);
                0.0     1.0  0.0;
            -sin(beta)  0.0  cos(beta)]

RotZ_g = [cos(gamma) -sin(gamma) 0.0;
            sin(gamma)  cos(gamma) 0.0;
            0.0     0.0    1.0] 

    prevCoordinates = fragment.Coords
    prevXAtom = fragment.XAtom
    pivotCoords = fragment.Coords[pivotAtomID,:]'
    prevCoordinates = prevCoordinates .- pivotCoords #This is meant to rotate around the origin
    prevXAtom = prevXAtom .- pivotCoords #This is meant to rotate around the origin

    newCoordinates = similar(prevCoordinates)
    newXAtom = similar(prevXAtom)
    nAtoms = size(prevCoordinates,1)
    #still missing rotation of the XAtom

    for kAtom in 1:nAtoms
        v = (prevCoordinates[kAtom,:]) #No need for transpose
        newCoordinates[kAtom,:] = transpose(RotZ_g * RotX_b * RotZ_a * v)
    end
    newXAtom = transpose(RotZ_g * RotX_b * RotZ_a * transpose(prevXAtom))

    newFragment = deepcopy(fragment)
    newFragment.Coords = newCoordinates .+ pivotCoords
    newFragment.XAtom = newXAtom .+ pivotCoords

    return newFragment
end

function attachFragment(growingFragment::Fragment,newFragment::Fragment)
    #This function takes two fragments (a growingFragment and a newFragment) and joins them, forming a bond using the definitions of ends of each fragment
    
    growingFragment.nAtomInternal += newFragment.nAtomInternal
    
    growingFragment.Masses = cat(growingFragment.Masses,newFragment.Masses,dims=1)
    
    currentCoords = newFragment.Coords
    currentCaps = newFragment.Caps

    #For now, beginAtom and endAtom are given manually, but they should be specified in the fragment information
    #and endAtom in growingFragment should be changing as the fragment grows.

    beginAtom = newFragment.idBegin        #9 - alcohol O
    endAtom = growingFragment.idEnd             #2 - carboxylic C
    nextEndAtom = growingFragment.idNextToEnd   #1 - carboxylic =O
    idCap = growingFragment.idCap               #1 - carboxylic O(H)

    #This bondVector is now done in an ad hoc way, but the way to find this vector should be based on geometries of bond angles and distances
    #I have to devise an algorithm to find this vector
    bondVector = growingFragment.Caps[1, :] - growingFragment.Coords[endAtom, :];
    
    #newTransformedCoords = ((newTransformedCoords .- transpose(newTransformedCoords[beginAtom, :])) 
    #                        .+ transpose(growingFragment.Coords[endAtom, :])) .+ transpose(bondVector)
    
    #This part is where the Rotation along angles and bonds should occur
    newTransformedCoords = ((currentCoords .- transpose(currentCoords[beginAtom, :])) 
                            .+ transpose(growingFragment.Caps[idCap, :]));

    # println("idCap is:")
    # println(idCap)
    # println(currentCaps)
    newTransformedCaps = ((currentCaps .- transpose(currentCoords[beginAtom, :])) .+ transpose(growingFragment.Caps[idCap, :]));
    # println(transpose(currentCoords[beginAtom, :]))
    # println(transpose(growingFragment.Caps[idCap, :]))
    # println(newTransformedCaps)

    growingFragment.Coords = cat(growingFragment.Coords,newTransformedCoords,dims=1);
    growingFragment.Caps = cat(growingFragment.Caps,newTransformedCaps,dims=1);

    
    growingFragment.idEnd += newFragment.nAtomInternal; #I have to subtract the previous Fragment (not the Growing Fragment) idEnd and add the newFragment idEnd, to account for different residues
    growingFragment.idNextToEnd += newFragment.nAtomInternal;
    growingFragment.idCap += newFragment.nAtomCaps
end

function attachSimpleFragment(growingFragment::SimpleFragment,newFragment::SimpleFragment)
    #This function takes two fragments (a growingFragment and a newFragment) and joins them, forming a bond using the definitions of ends of each fragment
    #It works with SimpleFragment:

    #fragment = SimpleFragment(nAtoms,        # nAtoms::Int
    #     idBegin,        # idBegin::Int
    #     idEnd, # idEnd::Int
    #     Masses, # Masses::Vector{Int64}
    #     Coords  # Coords::Matrix{Float64}
    # )
    
    #and therefore has less functionality
    #It doesn't have caps, which complicate the code, but it does need an X-atom, to form the new bond vector

    nAtoms = growingFragment.nAtoms

    growingFragment.Masses = cat(growingFragment.Masses,newFragment.Masses,dims=1)
    
    currentCoords = newFragment.Coords
    
    beginAtom = newFragment.idBegin        #9 - alcohol O
    endAtom = growingFragment.idEnd             #2 - carboxylic C

    #This bondVector is now done in an ad hoc way, but the way to find this vector should be based on geometries of bond angles and distances
    #I have to devise an algorithm to find this vector
    
    #XAtom are the coordinates of a H bonded to the idEnd atom, so the bondlenght should be scaled depending on the atom types
    bondVector = 1.4*(growingFragment.XAtom .- transpose(growingFragment.Coords[endAtom, :]));
    #println("bondVector is: ")
    #println(bondVector)

    #This part is where the Rotation along angles and bonds should occur
    newTransformedCoords = ((currentCoords .- transpose(currentCoords[beginAtom, :])) 
                            .+ transpose(growingFragment.Coords[endAtom, :])
                            .+ bondVector);
    newXAtomCoords = ((newFragment.XAtom .- transpose(currentCoords[beginAtom, :])) 
                            .+ transpose(growingFragment.Coords[endAtom, :])
                            .+ bondVector);

    growingFragment.Coords = cat(growingFragment.Coords,newTransformedCoords,dims=1);
 
    growingFragment.idEnd += newFragment.nAtoms; #I have to subtract the previous Fragment (not the Growing Fragment) idEnd and add the newFragment idEnd, to account for different residues
    growingFragment.XAtom = newXAtomCoords

    growingFragment.nAtoms = nAtoms + size(newFragment.Coords,1)

    return growingFragment
end



function rotateBondAngle(oldFragment::Fragment,pivotAtomInd,newAngle,topology)

    oldFragmentCoords = oldFragment.Coords
    rotatedFragment = oldFragment

    pivotRow = findall(x->x==pivotAtomInd,topology.angleIndList[:,2])
    prevPivotInd = Int64(first(topology.angleIndList[pivotRow,1]))
    nextPivotInd = Int64(first(topology.angleIndList[pivotRow,3]))
    currentAngle = topology.angleIndList[pivotRow,4]
    pivotAtomInd = Int64(pivotAtomInd)
    bond2_1 = oldFragment.Coords[prevPivotInd,:]-oldFragment.Coords[pivotAtomInd,:]
    bond2_3 = oldFragment.Coords[nextPivotInd,:]-oldFragment.Coords[pivotAtomInd,:]

    #This was wrong before, the vectorAxis must be a unit vector
    vectorAxis = cross(bond2_1,bond2_3)
    vectorAxis = vectorAxis/norm(vectorAxis)

    theta = (newAngle-first(currentAngle))*pi/180


    rodriguesRotations(rotatedFragment,pivotAtomInd,vectorAxis,theta)

    rotatedFragment.Coords[1:pivotAtomInd-1,:] = oldFragmentCoords[1:pivotAtomInd-1,:]

    return rotatedFragment

end
