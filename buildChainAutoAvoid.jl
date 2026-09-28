#buildChain.jl
#This script is based on the buildMolecule.jl script, but it is aimed to be more generally
#It uses fragments defined as XYZ coordinates

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
include("monomers.jl")

# #Definition of the monomer (I also have separate XYZ files for the starting and end monomer)
# Acrylamide_monomer = "C     0.000000000000      0.000000000000      0.000000000000
#                       H     0.000000000003      1.019689999998     -0.399256000006
#                       H    -0.883077999997     -0.509845000002     -0.399256000003
#                       C    -0.000000000011      0.000000000009      1.529988000000
#                       H     0.883077999986      0.509845000012      1.929244000003
#                       C    -0.000000000015     -1.433996037585      2.091464058406
#                       O    -0.000000000023     -1.719333811194      3.272742678927
#                       N    -0.000000000027     -2.541759345065      1.125151244354
#                       H    -0.378468417635     -2.312539998624      0.216113147949
#                       H    -0.378468417619     -3.411260018293      1.475665650290";

# #I need to define this atom as the missing place where the new fragment will get attached
# XAtom = [-0.898096987000      0.518516635000      1.831587572000]

# idBegin, idEnd = 1,4

# ################
# Acrlymide_end = "C     0.000000000000      0.000000000000      0.000000000000
#                 H     0.000000000000      1.019690000000     -0.399256000000
#                 H    -0.883078000000     -0.509845000000     -0.399256000000
#                 C     0.000000000000      0.000000000000      1.529988000000
#                 H     0.883078000000      0.509845000000      1.929244000000
#                 C     0.000000000000     -1.433996038000      2.091464058000
#                 O     0.000000000000     -1.719333811000      3.272742679000
#                 N     0.000000000000     -2.541759345000      1.125151244000
#                 H    -0.378468418000     -2.312539999000      0.216113148000
#                 H    -0.378468418000     -3.411260018000      1.475665650000
#                 H    -0.898096987000      0.518516635000      1.831587572000";

# AMP_monomer = "O     0.000000000000      0.000000000000      1.494856874000
# P     0.000000000000      0.000000000000      0.000000000000
# O    -1.157445167000     -0.818984107000     -0.420671871000
# O     1.330974569000     -0.489333563000     -0.426718464000
# O    -0.194549291000      1.505918781000     -0.489677166000
# C    -1.226937794000      2.313518213000      0.129482224000
# C    -1.813921733000      3.264959619000     -0.890358585000
# O    -0.922293080000      4.416590082000     -1.022907856000
# C    -1.952729638000      2.729863145000     -2.314166960000
# C    -0.587216502000      3.023914498000     -2.919928088000
# C    -0.219645207000      4.339729863000     -2.248888707000
# N     1.230128737000      4.478875027000     -1.953209699000
# C     2.107336568000      3.511113877000     -1.533231149000
# N     3.316110336000      3.949586076000     -1.354907866000
# C     3.236966918000      5.288894119000     -1.691718428000
# C     4.197577681000      6.323314632000     -1.709776747000
# N     5.473648493000      6.152278270000     -1.378890128000
# N     3.784111029000      7.542899392000     -2.085049982000
# C     2.507648993000      7.715270040000     -2.416194893000
# N     1.523029993000      6.843712932000     -2.434154772000
# C     1.967990022000      5.632454813000     -2.051968845000
# H    -2.018874412000      1.656420234000      0.516150238000
# H    -0.792441818000      2.890288894000      0.958200648000
# H    -2.800421479000      3.577782762000     -0.530227363000
# H    -0.526650301000      5.166283690000     -2.906166982000
# H     0.138664085000      2.233836185000     -2.676669906000
# H    -0.653341726000      3.136440025000     -4.013098218000
# H    -2.167336703000      1.649733165000     -2.310753863000
# H     1.817095091000      2.474298022000     -1.366143973000
# H     5.810952885000      5.223608980000     -1.088915122000
# H     6.126296775000      6.949268064000     -1.413580610000
# H     2.239895468000      8.721694671000     -2.715698896000";

# #I need to define this atom as the missing place where the new fragment will get attached
# #It should be the H atom attached to the endAtom
# XAtom = [0.447209972000     -0.778733265000      1.812051251000]
# #I change it slightly so that it is higher
# XAtom = [0.447209972000     -0.778733265000      2.012051251000]

# idBegin, idEnd = 9,1

################

idBegin = NIPAM_idBegin
idEnd = NIPAM_idEnd
XAtom = NIPAM_XAtom

atoms, coords, coordsXYZ = makeAtoms(NIPAM_monomer)                      
Masses = [ustrip(atoms[i].mass) for i in eachindex(atoms)]


GrowingMonomer = buildFragment(Masses,10*coordsXYZ,idBegin, idEnd, XAtom) #part of the program is configured in Angstrom
NewMonomer = buildFragment(Masses,10*coordsXYZ,idBegin, idEnd, XAtom)

#Rotation angles to rotate the new monomer along the z-axis, to avoid collision of amide fragments
alpha_0 = 120.0/180*pi
beta_0 = 0.0/180*pi
gamma_0 = 0.0/180*pi

#I need to generate a list of where each residue starts and ends in a Chain

Chain = []

#We initialize for the first residue
currResidueId = 1
currIdBegin = idBegin
currIdEnd = idEnd
currLastAtom = GrowingMonomer.nAtoms
push!(Chain,[currResidueId currIdBegin currIdEnd currLastAtom])

ntries = 500
alphastep = -2 #degrees
betastep = -2 #degrees
gammastep = -15 #degrees

    #these are initial values to test rotations
    global best_alpha = 0.0*pi/180.0
    global best_beta = 0.0
    global best_gamma = 0.0*pi/180.0
    kPivot = 1
    global best_pivotID = eligible_pivots[kPivot]
    worst_penalty = 10000

for k = 1:3
    
    #This attempts to find the rotation angles phi,theta, which are equivalent in my Euler matrices to alpha,beta
    #from the orientation of the vector joining the (current)idEnd and XAtom (which is supposed to change with changes in the GrowingMonomer)
    #the current idEnd must be calculated here, before adding the NewMonomer.nAtoms
    global XBond = transpose(GrowingMonomer.XAtom)-GrowingMonomer.Coords[currIdEnd,:]
    global nXBond = XBond/norm(XBond)
    global theta = acos(nXBond[3])
    global phi = acos(nXBond[1]/sin(theta))
    ##############################################

    #we first want to append the residueId, and new idBegin and idEnd
    global currResidueId = currResidueId + 1
    global currIdBegin = currIdBegin + NewMonomer.nAtoms
    global currIdEnd = currIdEnd + NewMonomer.nAtoms
    global currLastAtom = currLastAtom + NewMonomer.nAtoms

    push!(Chain,[currResidueId currIdBegin currIdEnd currLastAtom])

    global itry = 1
    #this loop will attempt to rotate around different bonds to avoid 
    
    global eligible_pivots = findall(x -> x == 12.0, NewMonomer.Masses)



    old_alpha = best_alpha
    old_beta = -(theta+30*pi/180) #This is because otherwise I have a very linear bond angle
    old_gamma = best_gamma
    old_pivotID = best_pivotID
    global old_penalty = worst_penalty

    while itry < ntries
        global new_alpha
        global new_beta
        global new_gamma
        global new_pivotID
        #This part is to iterate different values for rotations
        if mod(itry,3) == 1
            new_alpha = old_alpha + alphastep*pi/180
            new_beta = old_beta
            new_gamma = old_gamma
            new_pivotID = eligible_pivots[kPivot]
        elseif mod(itry,3) == 2
            new_alpha = old_alpha 
            new_beta = old_beta + (first(rand(1))-0.5)*betastep*pi/180
            new_gamma = old_gamma
            new_pivotID = eligible_pivots[kPivot]
        elseif mod(itry,3) == 0
            new_alpha = old_alpha 
            new_beta = old_beta 
            new_gamma = old_gamma + (first(rand(1))-0.5)*gammastep*pi/180
            new_pivotID = eligible_pivots[kPivot]
        # elseif mod(itry,4) == 3 #I want to have more chances of rotating angles before switching pivots
        #     new_alpha = old_alpha + alphastep*pi/180
        #     new_beta = old_beta
        #     new_pivotID = eligible_pivots[kPivot]
        end
        #This part does not work, because the Euler rotation rotates the whole fragment around pivotID,
        #not only the atoms including pivotID and following. So it rotates the fragment in an uncontrolled manner
        # if mod(itry,10) == 0 #Every 10 tries I change the pivot, and I reset the angles
        #     new_alpha = best_alpha
        #     new_beta = best_beta
        #     if kPivot < length(eligible_pivots)
        #         kPivot = kPivot + 1
        #     end
        #     new_pivotID = eligible_pivots[kPivot]
        # end
            # println("Current try number is: ")
            # println(itry)

            # println("Inside the while loop, these are the current values: ")
            # println("alpha, beta, gamma, pivotID ")
            # println([new_alpha, new_beta, new_gamma, new_pivotID])


        RotatedMonomer = eulerRotation(NewMonomer,new_pivotID,new_alpha,new_beta,new_gamma)
           

        global penalty = 0

        for k in eachindex(RotatedMonomer.Masses)
            local R1 = RotatedMonomer.Coords[k,:]
            for l in eachindex(GrowingMonomer.Masses)
                local R2 = GrowingMonomer.Coords[l,:]
                local dist12 = norm(R1-R2)
                if dist12 < 1.4 #this is in angstrom
                    penalty = penalty + 1
                end
                if dist12 < 1.0
                    penalty = penalty + 2
                end
                if dist12 < 0.5
                    penalty = penalty + 5
                end
            end
        end

        # println("Current penalty for close atoms is: ")
        # println(penalty)
        # println("Previous penalty for close atoms was: ")
        # println(old_penalty)


        if penalty < (old_penalty - 1.1)
            old_penalty = penalty
            old_alpha = new_alpha
            old_beta = new_beta
            old_gamma = new_gamma
            old_pivotID = new_pivotID
            #println("Current penalty is SMALLER than previous penalty")
            itry = itry + 1
        elseif penalty > old_penalty
            old_penalty = old_penalty
            old_alpha = old_alpha
            old_beta = old_beta
            old_gamma = old_gamma
            old_pivotID = old_pivotID
            #println("Current penalty is LARGER than previous penalty")
            itry = itry + 1
        elseif abs(penalty - old_penalty) < 1.1 #This should be to export out of the while loop the best values
            if (penalty < 10) || (itry == ntries - 1)
                global best_alpha = new_alpha
                global best_beta = new_beta
                global best_gamma = new_gamma
                global best_pivotID = new_pivotID
                #println("Current penalty is approximately equal to previous penalty, but it is also very small")
                break
            else
                old_penalty = penalty
                old_alpha = new_alpha
                old_beta = new_beta
                old_gamma = new_gamma
                old_pivotID = new_pivotID
                #println("Current penalty is approximately equal to previous penalty")
                itry = itry + 1
            end
        end

                global best_alpha = new_alpha
                global best_beta = new_beta
                global best_gamma = new_gamma
                global best_pivotID = new_pivotID
    end

    println("Fragment attachment number: ")
    println(k)
    println("After break of the while loop, these are the best values: ")
    println("itry, alpha, beta, gamma, pivotID ")
    println([itry, best_alpha, best_beta, best_gamma, new_pivotID])


    #RotatedMonomer = eulerRotation(NewMonomer,new_pivotID,new_alpha,new_beta,new_gamma)
    RotatedMonomer = eulerRotation(NewMonomer,best_pivotID,best_alpha,best_beta,new_gamma)
    global GrowingMonomer = attachSimpleFragment(GrowingMonomer,RotatedMonomer);
    
    # println("The last atom id in the growing fragment is")
    # println(growingFragment.idEnd)
    # println("growingFragment Coords are:")
    # print(growingFragment.Coords)
    # println("growingFragment Caps are:")
    # print(growingFragment.Caps)
end

bondIndList = findBondedAtoms(GrowingMonomer,elementTypeDict,atomRadiiDict,"angs")
#Based on this bondIndList and Chain information, I need to generate a script that finds 
#atoms which are in a residue m, are NOT currIdEnd for that residue, and are too close, or bound to 
#atoms in residue m-1, or even m-2
#after this, I need to make an algorithm to move them without creating further distortions.

Chain = reduce(vcat,Chain)

for resId in Chain[:,1]
    # println("resId: ")
    # println(resId)
    
    for atomId in Chain[resId,2]:Chain[resId,4]
        # println("atomId: ")
        # println(atomId)
    end
end

plotStructure(GrowingMonomer.Masses,GrowingMonomer.Coords,"angs");

#printXYZ(GrowingMonomer);

