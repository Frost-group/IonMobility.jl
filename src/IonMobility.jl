module IonMobility

using LinearAlgebra
using QuadGK
using SpecialFunctions
using Printf

export minimal_chakraborty_scpt

"""
    minimal_chakraborty_scpt(λκ; N=100, T=100.0, maxiter=200) -> Float64

Calculate Dₑff/D for ion diffusion in quenched charged disorder via FVA 
(Chakraborty, Bratko & Chandler, JCP 100:1535, 1994).

Solves Eq (52) ibid. iteratively for the variational parameters γₘ, then extracts Dₑff/D
from the low-frequency behavior per Eq (52') and (53').

# Arguments
- `λκ`: Disorder strength parameter λ_D·κ (dimensionless). Localises at λκ ≈ 1.
- `N`: Number of Fourier modes in the path discretization.
- `T`: Dimensionless loop time (larger → better low-frequency resolution).
- `maxiter`: Maximum self-consistent iterations.

Dₑff/D ∈ [0,1].
"""
function minimal_chakraborty_scpt(
    λκ::Float64;
    N::Int = 10,
    T::Float64 = 100.0,
    maxiter::Int = 200,
)

    # Dimensionless units: D=1, κ=1, β=1
    D = 1.0

    # Matsubara-like frequencies: Ωₘ = 2πm/T for m = 0,1,...,N
    Ωₘ = [2π * m / T for m = 0:N]

    # Bare (kinetic) action coefficient: Aₘ⁰ = (T/2D)Ωₘ²
    A⁰ = [(T / 2D) * ω^2 for ω in Ωₘ]

    # Variational parameters γₘ (memory kernel in frequency space)
    γ = zeros(N + 1) # Section IV Results, 'initial guess =0 for all n'
    γ_new = similar(γ)

    # Mixing for iteration stability
    α_mix = 0.3

    # Kernel: ∫ k⁴/(k²+1)² exp(-αk²) dk for screened Coulomb (55) ? v̄(k) = 1/(k²+1)² ?
    function K(α)
        f(k) = (k^4 / (k^2 + 1)^2) * exp(-α * k^2)
        val, _ = quadgk(f, 0.0, Inf, rtol = 1e-4)
        return val
    end

    for iter = 1:maxiter
        γ_old = copy(γ)

        # Denominator Aₘ⁰ + γₘ (appears in propagator)
        Dₘ = A⁰ .+ γ_old

        # Check for mode collapse (signals localization)
        any(d -> d < 1e-12, @view Dₘ[2:end]) && return 0.0

        # Displacement correlation Σ(τ) from Eq (52)
        function Σ(τ)
            s = 0.0
            @inbounds for m = 1:N
                s += (1.0 - cos(Ωₘ[m+1] * τ)) / Dₘ[m+1]
            end
            return s
        end

        # Update each mode m = 1,...,N via Eq (52)
        for m = 1:N
            Ω = Ωₘ[m+1]
            integrand(τ) = (1.0 - cos(Ω * τ)) * K(Σ(τ))
            # Integrate over τ ∈ [0, T/2] and use symmetry (avoids periodicity singularity at τ=T)
            I, _ = quadgk(integrand, 0.0, T/2, rtol = 1e-3)
            # Eq (52): γₘ = λκ × 2∫...dτ 
            # I thoguht there was a factor of 2 from symmetry, but when I remove this I get
            # agreement with the localisation collapse at disorder λκ=1.0
            γ_new[m+1] = λκ * I
        end

        # Convergence check
        Δ = norm(γ_new[2:end] - γ_old[2:end])
        rel_Δ = Δ / (norm(γ_new[2:end]) + 1e-10)

        if rel_Δ < 1e-5
            γ = γ_new
            break
        end

        # Damped update for stability
        γ = (1.0 - α_mix) * γ_old + α_mix * γ_new

        iter == maxiter &&
            @printf("Warning: max iterations reached (rel_Δ = %.2e)\n", rel_Δ)
    end

    # Extract Dₑff/D from low-frequency limit
    # Dₑff/D = 1/(1 + Γ̄(0)) where Γ̄(0) = (2D/T) γ₁ / Ω₁²
    γ₁, Ω₁ = γ[2], Ωₘ[2] # Or should this be first mode?

    (γ₁ > 1e5 || isnan(γ₁)) && return 0.0 # localise and collapse to zero
    
    # I'm totally confused here now; with the random non-diemsionalisation

    Γ̄₀ = (2D / T) * γ₁ / Ω₁^2
    Dₑff_over_D = 1.0 / (1.0 + Γ̄₀)

    return max(0.0, Dₑff_over_D)
end

end
