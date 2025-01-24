using LinearAlgebra
"""
Approximation de la solution du problème 

    min qₖ(s) = s'gₖ + 1/2 s' Hₖ s

        sous les contraintes s = -t gₖ, t > 0, ‖s‖ ≤ Δₖ

# Syntaxe

    s = cauchy(g, H, Δ; kwargs...)

# Entrées

    - g : (Vector{<:Real}) le vecteur gₖ
    - H : (Matrix{<:Real}) la matrice Hₖ
    - Δ : (Real) le scalaire Δₖ
    - kwargs  : les options sous formes d'arguments "keywords", c'est-à-dire des arguments nommés
        • tol_abs  : la tolérence absolue (optionnel, par défaut 1e-10)

# Sorties

    - s : (Vector{<:Real}) la solution du problème

# Exemple d'appel

    g = [0; 0]
    H = [7 0 ; 0 2]
    Δ = 1
    s = cauchy(g, H, Δ)

"""
function cauchy(g::Vector{<:Real}, H::Matrix{<:Real}, Δ::Real; tol_abs::Real = 1e-10)

    # Calculs des coefficients de la quadratique approchée
    a = (1/2) * g' * H * g
    b = -1 * g' * g 

    if norm(g, 2) == 0
        return 0
    elseif a == 0
        return (Δ / norm(g, 2))
    elseif a > 0
        # Calcul du minimum de la quadratique 
        t_min = -b / (2 * a)
    
        if t_min < 0
            t_min = 0
        elseif t_min > (Δ / norm(g, 2))
            t_min = (Δ / norm(g, 2))
        end
    elseif a < 0
        max = -b / (2 * a)

        if max >= (Δ / 2 * norm(g, 2)) 
            t_min = 0
        else
            t_min = (Δ / norm(g, 2))
        end
    end

    s = -1 * t_min * g
    return s
end
