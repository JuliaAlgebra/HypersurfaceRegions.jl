export membership


"""
    membership(chambers_list::Vector{Chamber}, p::Union{Vector{Float64},Vector{Int64}})

Input the list of chambers of the complement of hypersurface arrangements and a point p. 
Outputs the chamber the point p belongs to.

Options:
* `reltol = 1e-6`, `abstol = 1e-9`: parameters for the accuracy of the ODE solver.
```
"""
function membership(R::ChambersResult, p::T; kwargs...) where {T<:AbstractArray}
    f, f_denom, s = g(R)
    f_list = [f.expressions; f_denom]
    k = length(f_list)

    logg = sum(s[i] * log(f_list[i]^2) for i = 1:k)
    HC.differentiate(logg, variables(f))
    ∇logg = fixed(System(∇, variables = variables(f)))

    membership(R, p, ∇logg; kwargs...)
end

function membership(
    R::ChambersResult,
    p::Union{Vector{Float64},Vector{Int64}},
    ∇logg::AS;
    kwargs...,
) where {AS<:AbstractSystem}
    chambers_list = R.chamber_list

    i = _membership(R, p, ∇logg; kwargs...)

    if isnothing(i)
        return i
    else
        return chambers_list[i]
    end
end

function _membership(
    R::ChambersResult,
    p::Union{Vector{Float64},Vector{Int64}},
    ∇logg::AS;
    reltol::Float64 = 1e-6,
    abstol::Float64 = 1e-9,
    warning::Bool = true,
) where {AS<:AbstractSystem}
    f, _, _ = g(R)

    chambers_list = R.chamber_list
    critical_points_lists = [chamber.critical_points for chamber in chambers_list]
    critical_points = vcat(critical_points_lists...)

    evaluate_values = f(p)
    if any(abs.(evaluate_values) .< 1e-10)
        if warning
            println("The point is on one of the hypersurfaces")
        end
        return nothing
    end
    critical_point_index, failed_info =
        limit_critical_point(set_up_ode(∇logg), p, critical_points, reltol, abstol)
    if critical_point_index == -1
        return nothing
    end

    end_critical_point = critical_points[critical_point_index]
    for (i, critial_point_list) in enumerate(critical_points_lists)
        if end_critical_point in critial_point_list
            return i
        end
    end
end
