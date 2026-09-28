using Molly
using GLMakie
include("GeomUtils.jl")
import Random
#using .GeomUtils

#In this file I'm replicating mollyDiatomicTest.jl, but attempting to read an xyz file and use it to define atoms


fileString = "O     5.000000000000      5.000000000000      5.000000000000
              H     5.000000000008      5.759336999999      5.596043000001
              H     5.000000000008      4.240663000000      5.596042999999"

atoms,coords,coordsXYZ = makeAtoms(fileString)
#I think there is a problem with this, which is different from using place_atoms
n_atoms = size(atoms,1)
boundary = CubicBoundary(3.0u"nm")
temp = 298.0u"K"

coords = place_atoms_at_coordinates(n_atoms, coordsXYZ,boundary)

velocities = [random_velocity((atoms[i].mass), 0.2*temp) for i in 1:n_atoms]

#This should be taken to a specific makeBonds function in geomUtils.jl
n_bonds_tot = 2
bonds = InteractionList2Atoms(
    [1,1],           # First atom indices
    [2,3],           # Second atom indices
    [HarmonicBond(k=464.0u"kJ * mol^-1 * nm^-2", r0=0.1u"nm") for _ in 1:n_bonds_tot],
)

n_angles_tot = 1
angles = InteractionList3Atoms(
    2,
    1,
    3,
    [CosineAngle(k=30.0u"kJ * mol^-1", θ0=109.0*pi/180.0)]

)

specific_inter_lists = (bonds,angles,)

# All pairs apart from bonded pairs are eligible for non-bonded interactions
#This should also be done in a specific function that loops over bonded atoms (and also consider second and third neighbors)
eligible = trues(n_atoms, n_atoms)
for i in 1:n_atoms
    eligible[1, 2] = false
    eligible[2, 1] = false
    eligible[1, 3] = false
    eligible[3, 1] = false
end

neighbor_finder = DistanceNeighborFinder(
    eligible=eligible,
    n_steps=10,
    dist_cutoff=1.5u"nm",
)

cutoff = DistanceCutoff(1.2u"nm")
pairwise_inters = (LennardJones(use_neighbors=true, cutoff=cutoff),)

sys = System(
    atoms=atoms,
    coords=coords,
    boundary=boundary,
    velocities=velocities,
    pairwise_inters=pairwise_inters,
    specific_inter_lists=specific_inter_lists,
    neighbor_finder=neighbor_finder,
    loggers=(
        temp=TemperatureLogger(10),
        coords=CoordinatesLogger(10),
        energy=TotalEnergyLogger(10),
        #writer=StructureWriter(10, "sim_diatomic_b.xyz"),
    ),
)

simulator = VelocityVerlet(
    dt=0.001u"ps",
    coupling=AndersenThermostat(temp, 1.0u"ps"),
)
simulate!(sys, simulator, 10_000)


visualize(
    sys.loggers.coords,
    boundary,
    "sim_H2O_a.mp4";
    connections=[(1,2);(1,3)],
)