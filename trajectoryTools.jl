#This file contains functions to analyze trajectory files, as well as allxyz files from Orca
#Some of the functions allow calculation of Lennard-Jones and Coulomb energies from the atomic positions of each structure in the trajectory

using Chemfiles
using GLMakie
using LinearAlgebra
using Statistics


function readEnergies(filename::String)
        
        traj_file = open(filename)
        
        energies = []
        for ln in eachline(filename)
            #print(ln)
            if contains(ln,"energy") #this works for a trajectory file from Orca
                current_ln = ln
                current_fields = split(current_ln)
                println(current_fields)
                if contains(current_fields[1],"energy")
                    push!(energies,parse(Float64,current_fields[2]))
                end
            elseif contains(ln,"Relaxed Surface Scan Step") #this works for an allxyz file from a relaxed surface scan
                current_ln = ln
                current_fields = split(current_ln)
                println(current_fields)
                if contains(current_fields[10],"E")
                    push!(energies,parse(Float64,current_fields[11]))
                end
            end
        end
        return energies
end

function allxyz2xyz(filename::String)

    filenameFields = split(filename,".")
    if filenameFields[2] == "allxyz"
        newFilename = filenameFields[1] * ".xyz"
    end

    oldFileContent = open(filename);
    newFileContent = ""

        for ln in eachline(oldFileContent)
           
            if contains(ln,">") #this works for a trajectory file from Orca
                continue
            else
                newFileContent *=  ln * "\n"
            end
        end
    write(newFilename,newFileContent)
    #this generates a file with an extra empty line at the beginning which causes it not to be read by Trajectory
end

function readTrajectory(filename::String)

    copyBool = false

    filenameFields = split(filename,".")
    if filenameFields[2] == "allxyz"
        newFilename = filenameFields[1] * ".xyz"
        copyBool = true
    else
        newFilename = filename
        copyBool = false
    end

    #This copies the file to a new file with a new extension. It may be extremely inappropriate
    #if the trajectory files are large.
    if copyBool
        if isfile(newFilename) == false
            cp(filename, newFilename)

        end
    end

    trajectory = Trajectory(newFilename)
    n_steps = length(trajectory)
    frames = []
    atomPositions = []

        for k = 1:n_steps
            current_frame = read_step(trajectory,k-1)
            push!(frames,current_frame)
            #current_pos is a Julia array, and therefore starts counting from 1
            current_pos = positions(current_frame)
            push!(atomPositions,current_pos)
        end
    trajectory, frames, atomPositions

end

function readCharges(filename::String)
    #This function reads an Orca output file to retrieve Löwdin readCharges

    filenameFields = split(filename,".")
    if filenameFields[2] == "out"
        outFileContent = open(filename);
    else
        error("There is no Orca output file with this name.")
    end

    atomNames = []
    atomCharges = []
    atomSpins = []

    currentAtomNames = []
    currentAtomCharges = []
    currentAtomSpins = []

    accumValues = false
    pushValues = false

    for ln in eachline(outFileContent)
        
        if contains(ln,"LOEWDIN ATOMIC CHARGES AND SPIN POPULATIONS") #this works for an output file from Orca, as long as Loewdin charges have been printed for every geometry
            accumValues = true
            #I initialize again these arrays after each optimization done
            currentAtomNames = []
            currentAtomCharges = []
            currentAtomSpins = []
            continue
        end

        if contains(ln,"------")
            continue
        end

        if accumValues == true #contains(ln,":")
            lnContent = split(ln)
                if length(lnContent) == 5 #Let's make sure that in different calculation types this has always the same structure
                    push!(currentAtomNames,lnContent[2])
                    push!(currentAtomCharges,lnContent[4])
                    push!(currentAtomSpins,lnContent[5])
                end

            #println(lnContent)            

        end

        if isempty(ln) #an empty line, signaling the end of the Loewdin values
            accumValues = false
            pushValues = true
        end

            # while !contains(ln,"  ") #empty line signaling the end of the loewdin charges
            #     lnContent = split(ln)
            #     println(lnContent)
            # end
        
        if contains(ln,"*** OPTIMIZATION RUN DONE ***")

            push!(atomNames,currentAtomNames)
            push!(atomCharges,currentAtomCharges)
            push!(atomSpins,currentAtomSpins)

            currentAtomNames = []
            currentAtomCharges = []
            currentAtomSpins = []

            accumValues = false
            pushValues = false
        end
    end

    return     atomNames, atomCharges, atomSpins
end