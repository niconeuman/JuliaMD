using Molly
using GLMakie
GLMakie.activate!(px_per_unit = 2.0)
using Colors
import Random
include("generateFragment.jl")
include("manipulateFragment.jl")
include("plotStructure.jl")
include("topology.jl")

include("geomUtils.jl")
include("AtomParameters.jl")

#25sep2026
#In version 0.12, I want to implement the posibility to read the gaff2 forcefield parameters.
#For this the makeBonds, makeAngles, makeDihedrals functions should be able to read the output structures
#of readFFFile()



DiLacticAcid = "O         -12.47173996      -2.16439717      -3.31640626
                H         -12.24705194      -2.44379619      -4.22651933
                O         -11.36566187      -0.87383107      -5.46740242
                O         -10.62367582       0.85784507      -4.15404632
                C         -11.25663087      -0.32600003      -4.37520234
                C         -11.75310990      -0.95327507      -3.07083224
                H         -10.87651684      -1.21449709      -2.46918619
                C         -12.67021297      -0.01296500      -2.30820318
                H         -13.05780100      -0.49900704      -1.40651211
                H         -12.15434893       0.90613507      -2.01383315
                H         -13.53667104       0.26488202      -2.91899722
                O         -11.02790585       2.61555720      -7.31705456
                O         -12.53696896       2.16702817      -5.73760344
                C         -11.37548387       2.11794116      -6.11181147
                C         -10.19150378       1.53326612      -5.34516641
                H          -9.64959174       0.83239806      -5.99167946
                C          -9.25808571       2.65444520      -4.91073038
                H          -8.42172565       2.25370417      -4.32799633
                H          -8.85626668       3.20072325      -5.76975844
                H          -9.78097975       3.37052226      -4.26638833
                H         -11.87790491       2.91993122      -7.70094959"

TetraLacticAcid = "O   -12.471739960000     -2.164397170000     -3.316406260000
                H   -12.247051940000     -2.443796190000     -4.226519330000
                O   -11.365661870000     -0.873831070000     -5.467402420000
                O   -10.623675820000      0.857845070000     -4.154046320000
                C   -11.256630870000     -0.326000030000     -4.375202340000
                C   -11.753109900000     -0.953275070000     -3.070832240000
                H   -10.876516840000     -1.214497090000     -2.469186190000
                C   -12.670212970000     -0.012965000000     -2.308203180000
                H   -13.057801000000     -0.499007040000     -1.406512110000
                H   -12.154348930000      0.906135070000     -2.013833150000
                H   -13.536671040000      0.264882020000     -2.918997220000
                O   -11.027905850000      2.615557200000     -7.317054560000
                O   -12.536968960000      2.167028170000     -5.737603440000
                C   -11.375483870000      2.117941160000     -6.111811470000
                C   -10.191503780000      1.533266120000     -5.345166410000
                H    -9.649591740000      0.832398060000     -5.991679460000
                C    -9.258085710000      2.654445200000     -4.910730380000
                H    -8.421725650000      2.253704170000     -4.327996330000
                H    -8.856266680000      3.200723250000     -5.769758440000
                H    -9.780979750000      3.370522260000     -4.266388330000
                C   -10.160401309756      3.720567260394     -7.609652644552
                O    -8.977950741353      2.695213449194     -9.467187616001
                C    -9.086910126429      3.644112033034     -8.697364940962
                O    -8.242322774572      4.709463522917     -8.645118318325
                C    -7.250136631217      4.733813773223     -9.682763275754
                H    -9.661482244058      3.632863979714     -6.639076421426
                C   -10.954515147765      5.012560267941     -7.693327228295
                H   -11.751901176707      5.024698259736     -6.942636665318
                H   -10.318624695395      5.889840706991     -7.539910397590
                H   -11.436304554703      5.109718654471     -8.672846052875
                C    -7.870012318563      5.078670066913    -11.034882602653
                O    -7.001820615702      4.949929539506    -12.059925116301
                H    -7.537964930289      5.167357527790    -12.852258179807
                O    -9.017710917400      5.440877406206    -11.242446011151
                H    -6.758030789955      3.755757098868     -9.746800786592
                C    -6.228867862318      5.796096665503     -9.301172336080
                H    -5.801278374638      5.581355644051     -8.315919058577
                H    -5.415190280108      5.854257218233    -10.030797987334
                H    -6.698717335311      6.784198982669     -9.235930534157";

Phe = "C    -1.241946167370      4.419346038307      0.156876332166
        C     0.090478853034      5.063840981722      0.585443951737
        C     0.866544596340      5.518872556051     -0.665126191925
        N    -0.180944266750      6.227165010557      1.443710482997
        O     0.551282447039      6.549486228711     -1.237396117568
        O     1.971219165892      4.745252803077     -1.142072279074
        C    -2.018311606697      3.964163699088      1.407382886316
        H    -1.042122304748      3.564794371477     -0.473824172476
        H    -1.831805643302      5.141151625266     -0.390003467104
        H     0.680155975161      4.342349215277      1.132396607827
        H     0.692973523736      6.651158331813      1.723936510322
        H    -0.691189959533      5.929833669451      2.264229678814
        H     2.794179591476      5.170714190165     -0.889177306160
        C    -2.765382683411      4.891932576900      2.144007138562
        C    -3.471601180727      4.478213200117      3.280403954418
        C    -3.430903589967      3.136644907650      3.680191044338
        C    -2.684057177562      2.208808738929      2.943515280681
        C    -1.977882083105      2.622529101598      1.807102259925
        H    -2.797122615664      5.926798850518      1.835607461985
        H    -4.047815579343      5.193946061260      3.848593740209
        H    -3.975675501406      2.817499325387      4.556783244940
        H    -2.652768836070      1.173874456958      3.251802940133
        H    -1.402066806551      1.906582940804      1.238861390442";

boundary = CubicBoundary((3.0)u"nm")
temp = 298.0u"K"

ff, fields, AtomTypes, BondTypes, AngleTypes, TorsionTypes, NonBondedTypes = readFFFile("gaff2.xml")

#25sep2026: makeAtoms should be able to read the Lennard-Jones parameters
#However, the current implementation is based only on the element names
#If I want to use something like gaff2, I need to first assign atomTypes to all atoms in the structure
#For this I need to make the topology, and topology may or may not use the output of makeAtoms
#Either I make buildTopology work only on the xyz file, with no dependencies,
#or somehow I make a second function that assigns to the atoms, the correct Lennard-Jones parameters
#I also have a problem with charges, since these are specific for each molecule.
atoms, coords, coordsXYZ = makeAtoms(Phe)
n_atoms = length(atoms)

Masses = [ustrip(atoms[i].mass) for i in eachindex(atoms)]
velocities = [random_velocity((atoms[i].mass), 0.5*temp) for i in 1:n_atoms]
newFragment = buildFragment(Masses,coordsXYZ,1,n_atoms)

topology = buildTopology(newFragment,elementTypeDict,atomRadiiDict)
atomTypeArray = determineAtomTypes(atoms,topology)

bonds = makeBonds(atoms,topology,atomTypeArray,BondTypes)
angles = makeAngles(atoms,topology,atomTypeArray,AngleTypes)
torsions, torsionLog = makeDihedrals(atoms,topology,atomTypeArray,TorsionTypes)


specific_inter_lists = (bonds,angles,torsions,)



# All pairs apart from bonded pairs are eligible for non-bonded interactions
#This should also be done in a specific function that loops over bonded atoms (and also consider second and third neighbors)
eligible = trues(n_atoms, n_atoms)

#this is to remove all bonded pairs
for k in 1:size(topology.bondIndList,1)
        id1 = Int64(topology.bondIndList[k,1])
        id2 = Int64(topology.bondIndList[k,2])
        eligible[id1,id2] = false
        eligible[id2,id1] = false
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
    dt=0.0005u"ps",
    coupling=AndersenThermostat(temp, 1.0u"ps"),
)
simulate!(sys, simulator, 10_000)

n_bonds_tot = size(topology.bondIndList,1)

#fig = Figure(size = (1920, 1080))

#This is for plotting with specific colors
colors = []
markersizes = []
baseMarkerSize = 0.04
for l = 1:n_atoms
        element = keys(filter(p -> p.second[2] == Masses[l], elementTypeDict))
        current_color = (atomPlotParams[first(element)][2],atomPlotParams[first(element)][3],atomPlotParams[first(element)][4])
        current_markersize = atomPlotParams[first(element)][1]
        push!(colors,RGB(current_color[1],current_color[2],current_color[3]))
        push!(markersizes,baseMarkerSize*current_markersize)
end

# ax = Axis3(fig[1, 1], 
#     azimuth = π/4,       # Rotate 45 degrees
#     elevation = π/6,     # Elevate 30 degrees
#     viewmode = :fit)     # Prevent resizing during orientation changes

visualize(
    sys.loggers.coords,
    boundary,
    "sim_TetraLacticAcid.mp4";
    connections=[(Int64(topology.bondIndList[j,1]),Int64(topology.bondIndList[j,2])) for j in 1:n_bonds_tot],
    markersize = markersizes,
    color = colors,
    linewidth = 3.0,
)