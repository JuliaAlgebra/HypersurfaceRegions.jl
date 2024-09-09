export membership


"""
    membership(regions_list::Vector{Region}, p::Union{Vector{Float64},Vector{Int64}})

Input the list of regions of the complement of hypersurface arrangements and a point p. 
Outputs the region the point p belongs to.

Options:
* `reltol = 1e-6`, `abstol = 1e-9`: parameters for the accuracy of the ODE solver.
```
"""
function membership(R::RegionsResult, p::T; kwargs...) where {T<:AbstractArray}
    f, f_denom, s = g(R)
    f_list = [f.expressions; f_denom]
    k = length(f_list)

    logg = sum(s[i] * log(f_list[i]^2) for i = 1:k)
    ∇ = HC.differentiate(logg, HC.variables(f))
    ∇logg = fixed(System(∇, variables = HC.variables(f)))

    membership(R, p, ∇logg; kwargs...)
end

function membership(
    R::RegionsResult,
    p::Union{Vector{Float64},Vector{Int64}},
    ∇logg::AS;
    kwargs...,
) where {AS<:AbstractSystem}
  
    C = _membership(R, p, ∇logg; kwargs...)

    return C
end

function _membership(
    R::RegionsResult,
    p::Union{Vector{Float64},Vector{Int64}},
    ∇logg::AS;
    reltol::Float64 = 1e-6,
    abstol::Float64 = 1e-9,
    warning::Bool = true,
) where {AS<:AbstractSystem}
    f, _, _ = g(R)

    regions_list = R.region_list
    critical_points_lists = [region.critical_points for region in regions_list]
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
    for region in regions_list
        if end_critical_point in region.critical_points
            return region
        end
    end
end
