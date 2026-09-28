mutable struct Fragment
    nAtomInternal::Int
    nAtomCaps::Int
    idBegin::Int
    idEnd::Int
    idNextToEnd::Int
    Masses::Vector{Int64}
    Coords::Matrix{Float64}
    CapMasses::Vector{Int64}
    Caps::Matrix{Float64}
    idCap::Int
end

mutable struct SimpleFragment
    nAtoms::Int
    idBegin::Int
    idEnd::Int
    Masses::Vector{Int64}
    Coords::Matrix{Float64}
    XAtom::Matrix{Float64} #This is a new addition, because I need to know where to attach the following fragment
end

function buildMolecule(monomerName::String)

    LA = cmp(monomerName,"LacticAcid")
    ME = cmp(monomerName,"methylene")

    if LA == 0

        fragment = Fragment(
        9,
        3,
        1,
        9,
        1,
        [16;16;12;1;12;1;1;1;12],
        [0.000000000000      0.000000000000      0.000000000000
        0.448776238073     -1.606349461098      2.135644541491
        -0.164650998051      0.594354364163      1.301242843839
        0.625915478542      1.351083839476      1.376905475355
        -1.551165348943      1.210551800604      1.467072382851
        -1.722772952505      1.976545465533      0.705706622152
        -2.331105907083      0.450466401240      1.374779282526
        -1.659528780157      1.671679949657      2.452046139786
        0.129044696981     -0.457741294115      2.348236507262],
        [16;1;1],
        [-0.008364815213      0.038685359877      3.604276012748
        0.196775423734     -0.696491498364      4.206280614644
        -0.184415406560      0.665699588912     -0.675457357650],
        1
        )

    
    elseif ME == 0
        fragment = Fragment(
        3,
        2,
        1,
        1,
        1,
        [12;1;1],
        [0.000000000000      0.000000000000      0.000000000000
         0.422558506649     -0.938711163654     -0.363959877620
         0.601668461205      0.835301983155     -0.363959877673],
        [1;1],
        [0.000000000004      0.000000000032      1.091879632939
        -1.024226967857      0.103409180466     -0.363959877646],
        1
        )
    else
        error("This fragment is not yet defined")    
    end

    return fragment
end

function buildFragment(Masses::Vector{Float64},Coords::Matrix{Float64},idBegin::Int,idEnd::Int,XAtom::Matrix{Float64} = [0.0 0.0 0.0])

    nAtoms = size(Masses,1)
    

    fragment = SimpleFragment(nAtoms,        # nAtoms::Int
        idBegin,        # idBegin::Int
        idEnd, # idEnd::Int
        Masses, # Masses::Vector{Int64}
        Coords,  # Coords::Matrix{Float64}
        XAtom
    )

    return fragment
end
