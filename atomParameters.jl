
#This is a Dictionary of elements. The element symbol is the key and the value is a 
#vector containing the atomic number, atomic weight, and σ "nm" and ϵ "kJ * mol^-1" Lennard-Jones parameters
elementTypeDict = Dict(
    "H" => [1,1,0.250,0.065],
    "C" => [6,12,0.340,0.360],
    "N" => [7,14,0.325,0.710],
    "O" => [8,16,0.296,0.650],
    "P" => [15,31,0.3694,0.960228],

)

factor = 1.15
atomRadiiDict = Dict(
    #These covalent radii are given in pm, from Wikipedia. 
    #In case of C, the largest radius (76 pm for sp3) is given    
    "H" => [factor*31],
    "B" => [factor*84],
    "C" => [factor*76],
    "N" => [factor*71],
    "O" => [factor*66],   
    "F" => [factor*57],
    "Al" => [factor*121],
    "Si" => [factor*111],
    "P" => [factor*107],
    "S" => [factor*105],
    "Cl" => [factor*102],
)

atomPlotParams = Dict(
    #relativemarkersize, R, G, B
    "H" => [1,200/255,200/255,200/255],
    "C" => [2,100/255,100/255,100/255],
    "N" => [2,80/255,20/255,200/255],
    "O" => [2,200/255,20/255,80/255],

)
