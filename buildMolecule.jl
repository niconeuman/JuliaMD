using GLMakie

include("generateFragment.jl")
include("manipulateFragment.jl")
include("plotStructure.jl")
include("topology.jl")

Molecule = buildMolecule("methylene")
GrowingMolecule = buildMolecule("LacticAcid")

#println(Molecule.Coords)
#println("So far so good")

# plotStructure([Molecule.Masses;Molecule.CapMasses],[Molecule.Coords;Molecule.Caps])
# vectorAxis = [0.0;0.0;1.0]
# theta = 0.75*pi
# rodriguesRotations(Molecule,2,vectorAxis,theta)

#println(Molecule.Coords)

growingFragment = GrowingMolecule
newFragment = buildMolecule("LacticAcid")

for k = 1:3
    
    attachFragment(growingFragment,newFragment);
    
    # println("The last atom id in the growing fragment is")
    # println(growingFragment.idEnd)
    # println("growingFragment Coords are:")
    # print(growingFragment.Coords)
    # println("growingFragment Caps are:")
    # print(growingFragment.Caps)
end

plotStructure(growingFragment.Masses,growingFragment.Coords);


topology = buildTopology(growingFragment)

#Test to see if I can rotate a part of a molecule
pivotAtomInd = 10.0
newAngle = 145
#newAngle = first(currentAngle)

# oldFragmentCoords = growingFragment.Coords
# rotatedFragment = growingFragment

rotatedFragment = rotateBondAngle(growingFragment,pivotAtomInd,newAngle,topology)

# pivotRow = findall(x->x==pivotAtomInd,topology.angleIndList[:,2])
# prevPivotInd = Int64(first(topology.angleIndList[pivotRow,1]))
# nextPivotInd = Int64(first(topology.angleIndList[pivotRow,3]))
# currentAngle = topology.angleIndList[pivotRow,4]
# pivotAtomInd = Int64(pivotAtomInd)
# bond2_1 = growingFragment.Coords[prevPivotInd,:]-growingFragment.Coords[pivotAtomInd,:]
# bond2_3 = growingFragment.Coords[nextPivotInd,:]-growingFragment.Coords[pivotAtomInd,:]

# #This was wrong before, the vectorAxis must be a unit vector
# vectorAxis = cross(bond2_1,bond2_3)
# vectorAxis = vectorAxis/norm(vectorAxis)

# theta = (newAngle-first(currentAngle))*pi/180

# # #for now, this function rotates everything in the molecule around the pivot, not a part of the structure, as needed
# oldFragmentCoords = growingFragment.Coords
# rotatedFragment = growingFragment
# rodriguesRotations(rotatedFragment,pivotAtomInd,vectorAxis,theta)

# rotatedFragment.Coords[1:pivotAtomInd-1,:] = oldFragmentCoords[1:pivotAtomInd-1,:]

plotStructure(rotatedFragment.Masses,rotatedFragment.Coords);