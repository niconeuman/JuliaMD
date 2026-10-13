#This file will contain tools for reading forcefield files into an appropriate container
#It will also contain charge assignment functions


function assignCharges(atoms, coords, totalCharge, atomTypeArray, topology)

    #Atom is not mutable, that is why I cannot define the charge once the Atom has been initialized
    nAtoms = length(atoms)
    atomCharges = []
    chargeValue = 0.0
    #initialize atomic charges
        for k = 1:nAtoms
            
            if atomTypeArray[k]["element"] == "C"
                chargeValue = 0.0
            elseif atomTypeArray[k]["element"] == "H"
                if atomTypeArray[k]["class"] == "H-C"
                    chargeValue = -0.1
                elseif atomTypeArray[k]["class"] == "H-N"
                    chargeValue = 0.15
                elseif atomTypeArray[k]["class"] == "H-O"
                    chargeValue = 0.25
                else
                    chargeValue = 0.0
                end
            elseif atomTypeArray[k]["element"] == "O"
                if atomTypeArray[k]["class"] == "C=O"
                    chargeValue = -1.1
                elseif atomTypeArray[k]["class"] == "C-O-C"
                    chargeValue = -0.3
                elseif atomTypeArray[k]["class"] == "C-O-H"
                    chargeValue = -0.6
                else
                    chargeValue = 0.0
                end
                chargeValue = -1.0
            elseif atomTypeArray[k]["element"] == "N"
                chargeValue = -1.0
            end

            push!(atomCharges,chargeValue)
        end
        
        currentTotalCharge = sum(atomCharges)
        #println("currentTotalCharge is $currentTotalCharge")
        diffTotalCharge = totalCharge - currentTotalCharge
        #println("diffTotalCharge is $diffTotalCharge")
        atomCharges = atomCharges .+ diffTotalCharge/nAtoms
        currentTotalCharge = sum(atomCharges)
        #println("currentTotalCharge is $currentTotalCharge")

        #Charge equilibration algorithm


        #Initialize χᵢ and J\_ii
        
        χᵢ = zeros(nAtoms)
        Jᵢᵢ= zeros(nAtoms)

        for k = 1:nAtoms
            if atomTypeArray[k]["element"] == "C"
                χᵢ[k] = 5.343
                Jᵢᵢ[k] = 10.126
            elseif atomTypeArray[k]["element"] == "N"
                χᵢ[k] = 6.899
                Jᵢᵢ[k] = 11.760
            elseif atomTypeArray[k]["element"] == "O"
                χᵢ[k] = 8.741
                Jᵢᵢ[k] = 13.364
            elseif atomTypeArray[k]["element"] == "F"
                χᵢ[k] = 10.874
                Jᵢᵢ[k] = 14.948
            elseif atomTypeArray[k]["element"] == "S"
                χᵢ[k] = 6.928
                Jᵢᵢ[k] = 8.972
            elseif atomTypeArray[k]["element"] == "H"
                χᵢ[k] = 4.528
                Jᵢᵢ[k] = 13.8904
            else
                χᵢ[k] = 0.0
                Jᵢᵢ[k] = 0.0
                error("no QEq parameters for $(atomTypeArray[k]["element"])")               
            end                
            
        end

        #This will not change during the iteration
        χ₁ = χᵢ[1]

        #Initialize D vector
        D₀ = zeros(nAtoms)
        D₀[1] = -totalCharge
        for k = 2:nAtoms
            D₀[k] = (χᵢ[k] - χ₁)
        end

        #Initialize J Matrix
        J = zeros(nAtoms,nAtoms)
        damping = 0.1 #angstrom
        n = 2 #damping parameter
        kl = 0.4 #*0.529177 #Klondike parameter
        for k = 1:nAtoms
            for l = 1:nAtoms
                if k == l
                    J[k,l] = Jᵢᵢ[k]

                else
                    #In this algorithm R must be in Angstrom
                    Rₖₗ = ustrip(10*norm(coords[k,:]-coords[l,:]))

                    #limit when Rₖₗ → 0
                    Jₖₗ = 0.5*(Jᵢᵢ[k]+Jᵢᵢ[l])

                    #J[k,l] = 0*(Jₖₗ - 14.4/Rₖₗ)*exp(-(Rₖₗ/damping)^n) + 14.4/Rₖₗ

                    #Modification by DasGupta-Huzinaga
                    Rₖₗ = Rₖₗ/0.529177 #angstrom to bohr
                    #14.4 = 27.21*0.529177 is the conversion of Ha*Bohr to eV*A
                    #but because I have converted Rₖₗ to Bohr, I think I only need the 27.21 eV/Ha
                    J[k,l] = 27.21/(Rₖₗ+1/(0.5/27.21*(Jᵢᵢ[k]+Jᵢᵢ[l])*exp(kl*Rₖₗ)))

                end
            end
        end

        println(round.(J, digits = 2))        
        #Initialize C Matrix
        C = zeros(nAtoms,nAtoms)
        C[1,:] = ones(nAtoms)

        for k = 2:nAtoms
            for l = 1:nAtoms
                C[k,l] = J[k,l]-J[1,l]
            end
        end

        Q = C \ -D₀

        equilCharges = Q

        #Dₚᵣ = D₀
        #Dₙ = zeros(nAtoms)
        
        tolChange = 1e-2
        maxiter = 500
        deltaC = 0.01
        alpha_k = 0.5
        iter = 1
        #I need to solve for part of the C matrix (first row)
        #So I will use a steepest descent algorithm
        # while (norm(C*D₀+D₀) > tolChange && iter < maxiter)

        #     targetQuantity = norm(C*D₀+D₀)
        #     gradQuantity = zeros(nAtoms)
        #     currentQuantity = 0
        #     C_1_old = C[1,:]
            
        #     #calculation of gradient of target Quantity
        #     for l = 1:nAtoms


        #         C[1,l] = C_1_old[l]+deltaC

        #         currentQuantity = norm(C*D₀+D₀)

        #         gradQuantity[l] = (currentQuantity - targetQuantity)/(deltaC)
        #     end

        #     #Steepest descent direction
        #     dk = -gradQuantity

        #     #I have to update the charges using the gradient
        #     C_1_new = C_1_old + alpha_k*dk
        #     C[1,:] = C_1_new
        #     # if iter < 3
        #     #     println("currentQuantity is: ")
        #     #     println(currentQuantity)
        #     #     println("targetQuantity is: ")
        #     #     println(targetQuantity)
        #     #     println(C*D₀+D₀)
        #     # end
        #     iter += 1
        # end

        # equilCharges = C[1,:]
        # println(iter)

        return atomCharges,equilCharges

end

