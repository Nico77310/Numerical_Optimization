using LinearAlgebra
"""
Approximation de la solution du problème 

    min qₖ(s) = s'gₖ + 1/2 s' Hₖ s, sous la contrainte ‖s‖ ≤ Δₖ 


# Syntaxe

    s = gct(g, H, Δ; kwargs...)

# Entrées

    - g : (Vector{<:Real}) le vecteur gₖ, 
    - H : (Matrix{<:Real}) la matrice Hₖ
    - Δ : (Real) le scalaire Δₖ
    - kwargs  : les options sous formes d'arguments "keywords", c'est-à-dire des arguments nommés
        • max_iter : le nombre maximal d'iterations (optionnel, par défaut 100)
        • tol_abs  : la tolérence absolue (optionnel, par défaut 1e-10)
        • tol_rel  : la tolérence relative (optionnel, par défaut 1e-8)

# Sorties

    - s : (Vector{<:Real}) une approximation de la solution du problème

# Exemple d'appel

    g = [0; 0]
    H = [7 0 ; 0 2]
    Δ = 1
    s = gct(g, H, Δ)

"""
function q(g::Vector{<:Real}, H::Matrix{<:Real}, s::Vector{<:Real})
    return g' * s + (1/2) * s' * H * s
end 

function calc_σ(g::Vector{<:Real}, H::Matrix{<:Real}, s::Vector{<:Real}, p::Vector{<:Real}, Δ)
    a = p' * p
    b = 2 * p' * s
    c = (s' * s) - Δ^2
    Δ₂ = b^2 - (4 * a * c)

    if Δ₂ > 0 
        σ1 = (-b - sqrt(Δ₂))/(2 * a)
        σ2 = (-b + sqrt(Δ₂))/(2 * a)

        if q(g, H, s + σ1 * p) < q(g, H, s + σ2 * p)
            σ = σ1
        else 
            σ = σ2
        end
    
    elseif Δ₂ == 0
        σ = -b / (2 * a)
    
    else 
        σ = 0
    end 

    return σ
end

function gct(g::Vector{<:Real}, H::Matrix{<:Real}, Δ::Real;
    max_iter::Integer = 100,
    tol_abs::Real = 1e-10, 
    tol_rel::Real = 1e-8)

    j = 0
    g0 = g
    s = zeros(length(g))
    p = -g

    while ((j <= max_iter) && (norm(g, 2) > max(norm(g0, 2) * tol_rel, tol_abs)))
        κ = p' * H * p
        if κ <= 0 
            σ = calc_σ(g, H, s, p, Δ)
            return s + σ * p
        end 
        α = (g' * g) / κ
        if norm(s + (α * p), 2) >= Δ
            σ = calc_σ(g, H, s, p, Δ)
            return s + σ * p
        end 
        s += α * p
        g_prev = g
        g += α * H * p
        β = (g' * g)/(g_prev' * g_prev)
        p = -g + β * p
        j += 1
    end

    return s
end
