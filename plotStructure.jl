using Colors
include("AtomParameters.jl")

function plotStructure(moleculeMasses,moleculeCoords,units::String)

    if units == "nm"
        moleculeCoords = moleculeCoords.*10
    elseif units == "angs"
        moleculeCoords = moleculeCoords
    else
        error("Distance units not correctly specified, please use either 'nm' or 'angs'")
    end
    
    nAtoms = size(moleculeMasses,1)
    cdata = zeros(nAtoms,3)

    atomColorMap = [:lightgray,:green,:green,:green,:green,:darkgrey,:darkblue,:darkred,:green,:green,:green,:green,:green,:green,:orange] #this only works when the number of elements matches the maximum atomic number. If not there is some kind of scaling so that the highest number matches the last color
    moleculeNumbers = floor.(Int64,moleculeMasses/2 .+ 0.1) #This is a hack for finding the atom numbers from the masses. It will not work when many heavy elements are included
    moleculeNumbers[moleculeNumbers.==0].=1 #This fixes H atoms, which otherwise would have zero atomic number

    #Define scale of the plot
    margin = 20.0

    xlow = minimum(moleculeCoords[:,1])-margin
    xhigh = maximum(moleculeCoords[:,1])+margin
    ylow = minimum(moleculeCoords[:,2])-margin
    yhigh = maximum(moleculeCoords[:,2])+margin
    zlow = minimum(moleculeCoords[:,3])-margin
    zhigh = maximum(moleculeCoords[:,3])+margin

    #This code is no longer used
    for k in 1:nAtoms
        if moleculeMasses[k,1] == 1
            cdata[k,:] = [0.8 0.8 0.8];
        elseif moleculeMasses[k,1] == 12
            cdata[k,:] = [0.6 0.6 0.6];
        elseif moleculeMasses[k,1] == 16
            cdata[k,:] = [0.8 0.2 0.2];
        else
            cdata[k,:] = [0.2 0.8 0.2]; #generic green for undefined atoms
        end
    end

    #bondLengthThreshold = 1.6 #This should be later changed to account for atom types
    bond_radii = [0.35;0.5;0.6;0.74;0.88;0.84;0.76;0.71;0.63;0.62;1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0]

    # factor = 1.15
    # atomRadiiDict = Dict(
    # #These covalent radii are given in pm, from Wikipedia. 
    # #In case of C, the largest radius (76 pm for sp3) is given    
    # "H" => [factor*31],
    # "B" => [factor*84],
    # "C" => [factor*76],
    # "N" => [factor*71],
    # "O" => [factor*66],   
    # "F" => [factor*57],
    # "Al" => [factor*121],
    # "Si" => [factor*111],
    # "P" => [factor*107],
    # "S" => [factor*105],
    # "Cl" => [factor*102],
    # )
    
    makeBonds = zeros(nAtoms*(nAtoms-1),6)
    currentBondMarker = 1
    
    for k in 1:nAtoms
        atomNumber1 = moleculeNumbers[k]
        atom1 = moleculeCoords[k,:]
        for l in k+1:nAtoms
            atomNumber2 = moleculeNumbers[l]
            atom2 = moleculeCoords[l,:]
            bond1_2 = atom2-atom1
            bondLengthThreshold = (bond_radii[atomNumber1]+bond_radii[atomNumber2])
            if norm(bond1_2) < bondLengthThreshold
                makeBonds[currentBondMarker,:] = [atom1 1*(atom2-atom1)] #the scaling factor may not be necessary, it came from Octave's quiver3 
                currentBondMarker += 1
            end
        end
    end

    makeBonds = makeBonds[1:currentBondMarker,:]
    
        
    # println()
    # println(makeBonds)

    fig = Figure()
    ax = Axis3(fig[1, 1], xlabel = "x label", ylabel = "y label", title = "Title", aspect = :data, limits = (xlow, xhigh, ylow, yhigh, zlow, zhigh))
    arrows3d!(ax, makeBonds[:,1],makeBonds[:,2],makeBonds[:,3],makeBonds[:,4],makeBonds[:,5],makeBonds[:,6], 
            color = :gray, tiplength = 0.0, tailradius = 0.1)
    #scatter!(ax, moleculeCoords, markersize = 20*sqrt.(moleculeMasses .+ 1), color = moleculeNumbers, colormap = atomColorMap)
    meshscatter!(ax, moleculeCoords, markersize = 0.1*sqrt.(moleculeMasses .+ 1), color = moleculeNumbers, colormap = atomColorMap[1:maximum(moleculeNumbers)])
    
    # for k = 1:nAtoms
    #     text!(ax, moleculeCoords[k,1]-0.1,moleculeCoords[k,2]+0.1,moleculeCoords[k,3]+0.1,text = string(k), fontsize = 50)
    # end
    display(fig)
    
    return 1.0
    
end

function interactivePlot(frames,energies,atomPositions,atomPlotParams)

    #For the energies plot
    energies = (energies .- minimum(energies))*2625.5002   
    #mkrSizeEnergies = 10*ones(size(energies))

    #For the structures plot
    nStructs = length(atomPositions)
    nAtoms = size(atomPositions[1],2)

    #This is for plotting with specific colors
    colors = []
    markersizes = []
    baseMarkerSize = 0.5

    for l = 1:nAtoms
            element = name(frames[1][l-1]) #because Chemfiles objects start counting from 0
            current_color = (atomPlotParams[element][2],atomPlotParams[element][3],atomPlotParams[element][4])
            current_markersize = atomPlotParams[element][1]
            push!(colors,RGB(current_color[1],current_color[2],current_color[3]))
            push!(markersizes,baseMarkerSize*current_markersize)
    end

    fig = Figure(size = (1500, 800))

    
    #Define scale of the plot
    # margin = 5.0

    # xlow = minimum(moleculeCoords[:,1])-margin
    # xhigh = maximum(moleculeCoords[:,1])+margin
    # ylow = minimum(moleculeCoords[:,2])-margin
    # yhigh = maximum(moleculeCoords[:,2])+margin
    # zlow = minimum(moleculeCoords[:,3])-margin
    # zhigh = maximum(moleculeCoords[:,3])+margin

    ax = Axis3(fig[1, 1], title = "Structure", aspect = :data)
    ax2 = Axis(fig[1, 2], title = "Energy (kJ/mol)")

    nStruct_sl = Slider(fig[2, 1:2], range = 1:1:nStructs, startvalue = 1)
    idx_obs = nStruct_sl.value

    moleculeCoords = lift(nStruct_sl.value) do n
        atomPositions[n]
    end

    markerSizeEnergies = @lift(begin
        mkrSizeEnergies = fill(12, nStructs)
        mkrSizeEnergies[$idx_obs] = 20
        mkrSizeEnergies
    end)

    colorEnergies = @lift(begin
        colorsEn = fill(RGB(0.3,0.3,0.7), nStructs)
        colorsEn[$idx_obs] = RGB(0.8,0.8,0.4)
        colorsEn
    end)

    #println(markerSizeEnergies)
    #mkrSizeEnergies[highlightPoint[]] = 15

    meshscatter!(ax, moleculeCoords, markersize = markersizes, color = colors)
    scatter!(ax2,1:nStructs,energies, markersize = markerSizeEnergies, color = colorEnergies)
    lines!(ax2,1:nStructs,energies,linewidth = 2)
    fig

end