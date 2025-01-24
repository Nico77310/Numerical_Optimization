# Ecrire les tests de l'algorithme du pas de Cauchy
using Test

function tester_cauchy(cauchy::Function)

    #Tolérance utilisé dans les tests
    tol_test = 1e-3

	Test.@testset "Pas de Cauchy" begin
        #Test 1
        g = [1 ; 1]
        H = [3 1 ; 1 2]
        Δ = sqrt(2)
        s = cauchy(g, H, Δ)
        Test.@test s ≈ (-0.285714285714285)*g atol = tol_test

        #Test 2 
        g = [1 ; 1]
        H = [0 -1 ; -1 0]
        Δ = sqrt(2)
        s = cauchy(g, H, Δ)
        Test.@test s ≈ -1.0*g atol = tol_test
        
    end

end