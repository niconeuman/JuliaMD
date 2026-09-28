#This program reads trajectory (xyz or allxyz) and output files from Orca calculations in order to retrieve
#energies, structures, charges and other information needed to fit forcefield parameters against these trajectories
#It will also calculate contributions from Coulomb and Lennard-Jones terms, in order to achieve a better fit of torsions and angle force constants

include("geomUtils.jl")
include("trajectoryTools.jl")
include("plotStructure.jl")



#filename
allxyzFilename = "Biphenyl_pO_pO_Torsion_r2SCAN3c_haa.allxyz"
outputFilename = "Biphenyl_pO_pO_Torsion_r2SCAN3c_haa.out"
xyzFilename = "Biphenyl_pO_pO_Torsion_r2SCAN3c_haa.xyz"


energies = readEnergies(allxyzFilename)
allxyz2xyz(allxyzFilename)
trajectory, frames, atomPositions = readTrajectory(xyzFilename)
atomNames, atomCharges, atomSpins = readCharges(outputFilename)

#Functions to move to dedicate Julia file
sigma_ = []
epsilon_ = []
for k in 1:length(atomNames[1])
    push!(sigma_,elementTypeDict[atomNames[1][k]][3])
    push!(epsilon_,elementTypeDict[atomNames[1][k]][4])

end


function energyLJ(atomNames,sigma_,epsilon_,atomPositions)

    atomPositions .= atomPositions./10 #to transform from angstrom to nm

    if size(atomPositions,1) > 1 #This means if I have atomPositions for each step in a trajectory
        @assert length(sigma_) == size(atomPositions[1],2)
        nAtoms = length(sigma_)

        LJTrajList = [] #I want an array the size of the number of structures in the trajectory
        LJTrajEnergy = []
        for k = 1:size(atomPositions,1) #iterate over each structure
            
            LJList = []
            for i = 1:nAtoms
                for j = i+1:nAtoms
                    atomRadius_i = atomRadiiDict[atomNames[1][i]] #in pm
                    atomRadius_j = atomRadiiDict[atomNames[1][j]] #in pm
                    bondLengthThreshold_ij = (atomRadius_i + atomRadius_j)/1000 #pm to nm

                    sigma_ij = (sigma_[i] + sigma_[j])/2
                    epsilon_ij = sqrt(epsilon_[i]*epsilon_[j])
                    R1 = atomPositions[k][:,i]
                    R2 = atomPositions[k][:,j]
                    r12 = norm(R2-R1)

                    #This is a hack to avoid bonded pairs
                    if r12 > first(bondLengthThreshold_ij) #nm
                        r6 = (sigma_ij/r12)^6
                        r12 = r6^2                    
                        LJ_ij = 4*epsilon_ij*(r12-r6)
                        push!(LJList,[i  j  LJ_ij])
                    end
                end
            end
            push!(LJTrajList,LJList)
            push!(LJTrajEnergy,sum(reduce(vcat,LJList)[:,3]))
        end
    end

    return LJTrajList, LJTrajEnergy
end
