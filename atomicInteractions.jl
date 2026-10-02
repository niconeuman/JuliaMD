
#Energy calculations################################################################################
function calcCoulomb(frames,atomPositions,atomNames,atomCharges)
    TotalLJEnergy = 0.0

end

function nonBondedLJ(currentFragment::Fragment,topology)

    TotalLJEnergy = 0.0

    for k = 1:nAtom
        R1 = currentFragment.Coords[k,:]
        for l = k+1:nAtom
            R2 = currentFragment.Coords[l,:]
            LJPotential(Coord1,Coord2,σ,ϵ)
            TotalLJEnergy += LJPotential(R1,R2,σ,ϵ)
        end
    end

    return TotalLJEnergy

end

function LJPotential(Coord1,Coord2,σ,ϵ)

    r12 = norm(Coord2-Coord1)
    
    rr12 = σ/r12

    rr12_6 = rr12^6
    rr12_12 = rr12_6^2



    LJEnergy = 4*ϵ*(rr12_12-rr12_6)

    return LJEnergy
end

