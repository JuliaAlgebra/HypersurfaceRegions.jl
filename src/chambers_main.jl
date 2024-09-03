export chambers


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
# homogenize f using x_0 and then set x_0=0 x_1=1 to get an affine polynomial at infinity
function get_f_infty(
    f::Expression,
    var_list::Vector{Variable},
    selected_variable::Variable = var_list[1],
)
    exponents, coeffs = exponents_coefficients(f, var_list)
    total_degf_each_term = sum(exponents, dims = 1)
    degf = maximum(total_degf_each_term)
    poly = Expression(0)
    num_monomials = length(total_degf_each_term)
    for i = 1:num_monomials
        if total_degf_each_term[i] == degf
            exponent = exponents[:, i]
            coeff = coeffs[i]
            monomial = coeff
            for j = 1:length(exponent)
                monomial = monomial * var_list[j]^exponent[j]
            end
            poly += monomial
        end
    end

    result = get_f_affine(poly, var_list, selected_variable)
    if HC.degree(result) == 0
        # need to take log later so we need to make sure the output is not 1 or -1
        result = Expression(2)
    end

    return result
end



# input a point a at infinity (P^n-R^n)
# output a point in the affine linear space in some unbounded chamber of R^n-V(f)
function point_unbounded(f::Expression, a::Array{T}) where {T<:Real}
    new_a = vcat(1, a)
    @unique_var t
    f_t = subs(f, HC.variables(f) => t * new_a)
    S = HC.solve([f_t], t; show_progress = false)
    R = real_solutions(S)

    if isempty(R)
        t = 1.0
    else
        m = maximum([abs(maximum(R)[1]), abs(minimum(R)[1])])
        t = 5.0 * (m + 1) # relative increase of m
    end
    return t * new_a
end



"""
    chambers(f::Vector{Expression})
    chambers(f::System)

Input a list of hypersurfaces 'f = [f_1,...f_k]'.
Outputs the chambers in the complement of the hypersurface arrangement, whether they are bounded or not, their sign patterns, Euler characteristic and the indices of the critical points in each chamber.

Options:
* `show_progress = true`: if true, prints the progress of the computation to the terminal.
* `projective_fusion = true`: if `true`, the algorithm computes which of the chambers are fused at infinity.
* `s`: a list of integers `[s_1, ..., s_k, s_{k+1}]` such that `s_1, ..., s_k>0, s_{k+1}<0` and `2 s_{k+1} > s_1 deg(f_1) + ... + s_k deg(f_k)`.
* `epsilon = 1e-6`: how close from each critical point do we do the path tracking.
* `reltol = 1e-6`, `abstol = 1e-9`: parameters for the accuracy of the ODE solver.
* `monodromy_options = MonodromyOptions(max_loops_no_progress = 25)`: pass options for [monodromy](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/monodromy/).
* `start_pair_using_newton::Bool = false`: if true, the algorithm tries to compute a start pair for monodromy by using Newton's methods. Can reduce the number of critical points, but is less stable.


##  Example
```julia
using Chambers
@var x y
f = [x^2 + y^2 - 1; x^2 + y^2 - 4];
chambers(f)
```
"""
chambers(f::Vector{Expression}; kwargs...) = chambers(System(f); kwargs...)
function chambers(f::System; show_progress::Bool = true, kwargs...)


    #  ProgressMeter
    if show_progress
        progress = ChambersProgress(
            PM.ProgressUnknown(
                dt = 0.75,
                desc = "Computing chambers...",
                enabled = true,
                spinner = true,
            ),
        )
    else
        progress = nothing
    end

    sleep(1.0)
    update_progress!(progress)
    R = _chambers(f, progress; kwargs...)
    finish_progress!(progress)

    R
end


function _chambers(
    f::System,
    progress::Union{Nothing,ChambersProgress};
    projective_fusion::Bool = true,
    s::Union{Nothing,Vector{T}} = nothing,
    epsilon::Float64 = 1e-6,
    reltol::Float64 = 1e-6,
    abstol::Float64 = 1e-9,
    monodromy_options = HC.MonodromyOptions(max_loops_no_progress = 10),
    start_pair_using_newton::Bool = false,
    seed = nothing,
    kwargs...,
) where {T<:Real}


    # Stage 1: computing affine chambers
    set_stage!(progress, 1)
    #

    affine_output = _affine_chambers(
        f,
        progress;
        s = s,
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

    # Stage 2: computing projective chambers
    set_stage!(progress, 2)
    #

    # get polynomials at the infinity space
    poly_list_infty = []
    poly_list = f.expressions
    variable_list = f.variables
    poly_list_infty = map(fᵢ -> get_f_infty(fᵢ, variable_list), poly_list)

    if all(HC.degree.(poly_list_infty) .== 0)
        critical_points_infty = Vector{Vector{ComplexF64}}()
    else
        infty_output = _affine_chambers(
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
            return nothing
        end

        critical_points_infty = map(infty_output.chamber_list) do chamber
            chamber.critical_points[1]
        end

    end

    # Stage 3: connecting critical points at infinity to affine chambers
    set_ncritical_points!(progress, 0)
    set_stage!(progress, 3)
    #

    # get one point from each chamber
    R = affine_output.chamber_list
    N = nchambers(affine_output)

    graph = LG.SimpleGraph(N)

    prod_f = prod(poly_list)

    f, f_denom, s = g(affine_output)
    f_list = [f.expressions; f_denom]
    k = length(f_list)

    logg = sum(s[i] * log(f_list[i]^2) for i = 1:k)
    ∇ = HC.differentiate(logg, variable_list)
    ∇logg = fixed(System(∇, variables = variable_list))

    set_ncritical_points!(progress, length(critical_points_infty))
    j = 0

    unbounded = Vector{Int}()
    for critical_point in critical_points_infty

        unbounded_point = point_unbounded(prod_f, critical_point)
        C1 = _membership(affine_output, unbounded_point, ∇logg)
        chamber_1_index = number(C1)

        if !isnothing(chamber_1_index)
            C2 = _membership(affine_output, -unbounded_point, ∇logg)
            chamber_2_index = number(C2)

            if !isnothing(chamber_2_index)
                LG.add_edge!(graph, chamber_1_index, chamber_2_index)
                append!(unbounded, chamber_1_index)
                append!(unbounded, chamber_2_index)
            end
        end

        j += 1
        set_ncritical_points_classified!(progress, j)

        #@warn "Some critical points could not be associated to one chamber."
    end

    bounded = setdiff(1:N, unbounded)

    unique!(bounded)
    unique!(unbounded)

    R_new = Vector{Chamber}()
    for i in bounded
        Rᵢ = R[i]
        push!(R_new, Chamber(Rᵢ.sign, Rᵢ.χ, Rᵢ.μ, Rᵢ.critical_points, Rᵢ.g, true, Rᵢ.chamber_number))
    end
    for i in unbounded
        Rᵢ = R[i]
        push!(R_new, Chamber(Rᵢ.sign, Rᵢ.χ, Rᵢ.μ, Rᵢ.critical_points, Rᵢ.g, false, Rᵢ.chamber_number))
    end

    if projective_fusion

        partition = LG.connected_components(graph)
        projective_chambers = []
        for par in partition
            if in(par[1], unbounded)
                push!(projective_chambers, par)
            end
        end

        for b in bounded
            push!(projective_chambers, [b])
        end
    else
        projective_chambers = nothing
    end

    return ChambersResult(
        R_new,
        nchambers(affine_output),
        ncritical_complex(affine_output),
        ncritical_real(affine_output),
        projective_chambers,
        g(affine_output),
    )




end
