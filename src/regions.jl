export regions


# input a homogeneous polynomial f with variables x_0...x_n 
# homogenize f using x_0=1 
function get_f_affine(
    f::Expression,
    var_list::Vector{Variable},
    selected_variable::Variable = var_list[1],
)

    return subs(f, selected_variable => 1.0)
end


# input a polynomial f with variables x_1 ...x_n 
# homogenize f using x_0 and then set x_1=1 to get an affine polynomial at infinity
function get_f_infty(
    f::Expression,
    var_list::Vector{Variable},
    x0::Variable,
    selected_variable::Variable = var_list[1],
)
    exponents, coeffs = exponents_coefficients(f, var_list)
    total_degf_each_term = sum(exponents, dims = 1)
    degf = maximum(total_degf_each_term)

    poly = Expression(0)
    num_monomials = length(total_degf_each_term)

    for i = 1:num_monomials
        exponent = exponents[:, i]
        coeff = coeffs[i]
        monomial = coeff
        for j = 1:length(exponent)
            monomial = monomial * var_list[j]^exponent[j]
        end

        k = degf - total_degf_each_term[i]
        poly += monomial * x0^k
    end



    result = get_f_affine(poly, var_list, selected_variable)

    if HC.degree(result) == 0
        # need to take log later so we need to make sure the output is not 1 or -1
        result = Expression(2)
    end

    return result
end


function get_points_at_infinity(
    poly_list_infty,
    variable_list,
    progress,
    s,
    epsilon,
    reltol,
    abstol,
    monodromy_options,
    start_pair_using_newton,
    seed,
)
    #
    if all(HC.degree.(poly_list_infty) .== 0)
        critical_points_infty = [randn(Float64, length(variable_list) - 1)]
    else
        infty_output = _affine_regions(
            System(poly_list_infty, variables = variable_list[2:end]),
            progress;
            s = s,
            epsilon = epsilon,
            reltol = reltol,
            abstol = abstol,
            monodromy_options = monodromy_options,
            start_pair_using_newton = start_pair_using_newton,
            seed = seed,
        )

        if isnothing(infty_output)
            critical_points_infty = [randn(Float64, length(variable_list) - 1)]
        end

        critical_points_infty = map(infty_output.region_list) do region
            region.critical_points[1]
        end

    end
    critical_points_infty

end


# input a point a at infinity (P^n-R^n)
# output a point in the affine linear space in some unbounded region of R^n-V(f)
function point_unbounded(f::Expression, a::Array{T}, δ) where {T<:Real}
    new_a = vcat(1, a)
    λ = norm(new_a)
    new_a_normed = new_a ./ λ

    @unique_var t
    f_t = subs(f, HC.variables(f) => t * new_a_normed) |> expand
    # e, c = exponents_coefficients(f_t, [t])
    # LA.normalize!(c)

    # f_t = sum(cᵢ * t^eᵢ for (eᵢ, cᵢ) in zip(e, c) if abs(cᵢ) > 1e-15) # remove almost zero terms

    if degree(f_t) == 0
        t = 1.0
    else
        S = HC.solve([f_t], t; show_progress = false)
        R = first.(real_solutions(S))

        invδ = inv(δ)
        # only work outside the strip around infinity
        if δ == 0.0
            filter!(r -> abs(r / λ) < invδ, R)
        elseif δ > 0
            filter!(r -> r / λ < invδ && r > 0, R)
        elseif δ < 0
            filter!(r -> r / λ > invδ && r < 0, R)
        end

        if isempty(R)
            t = 1.0
        else

            if δ == 0.0
                m = max(abs(maximum(R)), abs(minimum(R)))
                t = 5.0 * (m + 1) # relative increase of m
            elseif δ > 0
                m = maximum(R)
                t = 5.0 * (m + 1) # relative increase of m
                if t / λ > invδ
                    t = sqrt(invδ * m) * λ
                end
            elseif δ < 0
                m = minimum(R)
                t = 5.0 * (m - 1) # relative decrease of m
                if t / λ < invδ
                    t = -sqrt(abs(invδ) * abs(m)) * λ
                end
            end

        end
    end

    return t * new_a_normed
end



"""
    regions(f::Vector{Expression})
    regions(f::System)

Input a list of hypersurfaces 'f = [f_1,...f_k]'.
Outputs the regions in the complement of the hypersurface arrangement, whether they are bounded or not, their sign patterns, Euler characteristic and the indices of the critical points in each region.

Options:
* `δ::Float64 = 1e-12`: Parameter that defines the strip around infinity.
* `target_parameters`: Specify parameters of the [System](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/systems/) `f` (if its has any).
* `show_progress = true`: if true, prints the progress of the computation to the terminal.
* `projective_fusion = true`: if `true`, the algorithm computes which of the regions are fused at infinity.
* `s`: exponents of the Morse function `f_1^(s_1) * ... * f_k^(s_k) * q^(s_k+1)`. Here, `s` is a list of integers `[s_1, ..., s_k, s_{k+1}]` such that `s_1, ..., s_k>0, s_{k+1}<0` and `2 s_{k+1} > s_1 deg(f_1) + ... + s_k deg(f_k)`.
* `epsilon = 1e-6`: how close from each critical point do we do the path tracking.
* `reltol = 1e-6`, `abstol = 1e-9`: parameters for the accuracy of the ODE solver.
* `monodromy_options = MonodromyOptions(max_loops_no_progress = 25)`: pass options for [monodromy](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/monodromy/).
* `start_pair_using_newton::Bool = false`: if true, the algorithm tries to compute a start pair for monodromy by using Newton's methods. Can reduce the number of critical points, but is less stable.


##  Example
```julia
using HypersurfaceRegions
@var x y
f = [x^2 + y^2 - 1; x^2 + y^2 - 4];
regions(f)
```

## Example with options 
```julia
regions(f; δ = 1e-12, 
            monodromy_options = MonodromyOptions(max_loops_no_progress = 20))
```
"""
regions(f::Vector{Expression}; kwargs...) = regions(System(f); kwargs...)
function regions(f::System; show_progress::Bool = true, kwargs...)


    #  ProgressMeter
    if show_progress
        progress = RegionsProgress(
            PM.ProgressUnknown(
                dt = 0.75,
                desc = "Computing regions...",
                enabled = true,
                spinner = true,
            ),
        )
    else
        progress = nothing
    end

    sleep(1.0)
    update_progress!(progress)
    R = _regions(f, progress; kwargs...)
    finish_progress!(progress)

    R
end


function _regions(
    f0::System,
    progress::Union{Nothing,RegionsProgress};
    δ::Float64 = 1e-12,
    target_parameters::Union{Nothing,Vector{T1}} = nothing,
    s::Union{Nothing,Vector{T}} = nothing,
    epsilon::Float64 = 1e-6,
    reltol::Float64 = 1e-6,
    abstol::Float64 = 1e-9,
    monodromy_options = HC.MonodromyOptions(max_loops_no_progress = 10),
    start_pair_using_newton::Bool = false,
    projective_fusion::Bool = true,
    seed = nothing,
    kwargs...,
) where {T<:Real,T1<:Number}

    if isnothing(f0)
        f = f0
    else
        f = System(f0(HC.variables(f0), target_parameters), variables = HC.variables(f0))
    end
    s = set_up_s(s, f)

    ####
    # Stage 1: computing affine regions
    set_stage!(progress, 1)
    #

    affine_output = _affine_regions(
        f,
        progress;
        s = s,
        target_parameters = nothing,
        epsilon = epsilon,
        reltol = reltol,
        abstol = abstol,
        monodromy_options = monodromy_options,
        start_pair_using_newton = start_pair_using_newton,
        seed = seed,
    )

    if isnothing(affine_output)
        return nothing
    end

    ####
    # Stage 2: compute critical points at infinity
    set_stage!(progress, 2)
    #

    # get polynomials at the infinity space
    poly_list = f.expressions
    variable_list = f.variables
    @unique_var x0
    f0 = map(fᵢ -> get_f_infty(fᵢ, variable_list, x0), poly_list)
    F0 = System(f0, variables = variable_list[2:end], parameters = [x0])

    cpt, _ = compute_critical_points(
        F0,
        s,
        monodromy_options,
        progress,
        start_pair_using_newton;
        target_parameters = [[0.0], [δ], [-δ]],
        kwargs...,
    )

    # critical points at infinity
    critical_points_infty = real_solutions(first(cpt[1]))
    # critical points at one side of the strip
    critical_points_infty_1 = real_solutions(first(cpt[2]))
    # critical points at the other side of the strip
    critical_points_infty_2 = real_solutions(first(cpt[3]))

    ####
    # Stage 3: connecting critical points at infinity to affine regions
    set_ncritical_points!(progress, 0)
    set_stage!(progress, 3)
    #

    # get one point from each region
    R = affine_output.region_list
    N = nregions(affine_output)

    graph = LG.SimpleGraph(N)

    prod_f = prod(poly_list)

    f, f_denom, s = g(affine_output)
    f_list = [f.expressions; f_denom]
    k = length(f_list)

    logg = sum(s[i] * log(f_list[i]^2) for i = 1:k)
    ∇ = HC.differentiate(logg, variable_list)
    ∇logg = fixed(System(∇, variables = variable_list))

    nall_crit =
        length(critical_points_infty) +
        length(critical_points_infty_1) +
        length(critical_points_infty_2)
    set_ncritical_points!(progress, nall_crit)
    j = 0

    unbounded = Vector{Int}()
    undecided = Vector{Int}()
    bounded = Vector{Int}()

    # critical points at infinity
    for critical_point in critical_points_infty
        unbounded_point = point_unbounded(prod_f, critical_point, 0.0)
        C1 = _membership(affine_output, unbounded_point, ∇logg)

        if !isnothing(C1)
            region_1_index = number(C1)
            append!(unbounded, region_1_index)

            if projective_fusion
                C2 = _membership(affine_output, -unbounded_point, ∇logg)

                if !isnothing(C2)
                    region_2_index = number(C2)

                    LG.add_edge!(graph, region_1_index, region_2_index)
                    append!(unbounded, region_2_index)
                end
            end
        end

        j += 1
        set_ncritical_points_classified!(progress, j)
    end

    # critical points at the strip around infinity
    for critical_point in critical_points_infty_1
        unbounded_point = point_unbounded(prod_f, critical_point, δ)
        C = _membership(affine_output, unbounded_point, ∇logg)

        if !isnothing(C)
            if !in(number(C), unbounded)
                append!(undecided, number(C))
            end
        end

        j += 1
        set_ncritical_points_classified!(progress, j)
    end
    for critical_point in critical_points_infty_2
        unbounded_point = point_unbounded(prod_f, critical_point, -δ)
        C = _membership(affine_output, unbounded_point, ∇logg)

        if !isnothing(C)
            if !in(number(C), unbounded)
                append!(undecided, number(C))
            end
        end

        j += 1
        set_ncritical_points_classified!(progress, j)
    end


    still_to_decide = setdiff(collect(1:N), vcat(unbounded, undecided))
    for i in still_to_decide
        C = R[i]
        c = critical_points(C) |> first
        invc = inv(c[1])

        if invc > δ || invc < -δ
            append!(bounded, number(C))
        else
            append!(undecided, number(C))
        end
    end

    unique!(bounded)
    unique!(unbounded)
    unique!(undecided)

    R_new = Vector{Region}()
    for i in unbounded
        Rᵢ = R[i]
        push!(
            R_new,
            Region(Rᵢ.sign, Rᵢ.χ, Rᵢ.μ, Rᵢ.critical_points, Rᵢ.g, 0, Rᵢ.region_number),
        )
    end
    for i in bounded
        Rᵢ = R[i]
        push!(
            R_new,
            Region(Rᵢ.sign, Rᵢ.χ, Rᵢ.μ, Rᵢ.critical_points, Rᵢ.g, 1, Rᵢ.region_number),
        )
    end
    for i in undecided
        Rᵢ = R[i]
        push!(
            R_new,
            Region(Rᵢ.sign, Rᵢ.χ, Rᵢ.μ, Rᵢ.critical_points, Rᵢ.g, 2, Rᵢ.region_number),
        )
    end


    if projective_fusion
        projective_regions = LG.connected_components(graph)
    else
        projective_regions = nothing
    end

    return RegionsResult(
        R_new,
        nregions(affine_output),
        ncritical_complex(affine_output),
        ncritical_real(affine_output),
        projective_regions,
        g(affine_output),
    )




end
