export affine_chambers


"""
    affine_chambers(f::Vector{Expression})
    affine_chambers(f::System)

Input a list of affine hypersurfaces 'f = [f_1,...f_k]'. 
Outputs the chambers in the complement of the hypersurface arrangement, their sign patterns, Euler characteristic and the indices of the critical points in each chamber.
Accepts the same options as `chambers`.

##  Example
```julia
using Chambers
@var x y
f = [x^2 + y^2 - 1; x^2 + y^2 - 4];
affine_chambers(f)
```
"""
affine_chambers(f::Vector{Expression}; kwargs...) = affine_chambers(System(f); kwargs...)
function affine_chambers(f::System; show_progress::Bool = true, kwargs...)

    #  ProgressMeter
    if show_progress
        progress = ChambersProgress(
            PM.ProgressUnknown(
                dt = 1.0,
                desc = "Computing chambers...",
                enabled = true,
                spinner = true,
            ),
        )
    else
        progress = nothing
    end

    set_stage!(progress, 1)
    update_progress!(progress)
    R = _affine_chambers(f, progress; kwargs...)
    finish_progress!(progress)

    R
end


function _affine_chambers(
    f::System,
    progress::Union{Nothing,ChambersProgress};
    s::Union{Nothing,Vector{T}} = nothing,
    epsilon::Float64 = 1e-6,
    reltol::Float64 = 1e-6,
    abstol::Float64 = 1e-9,
    monodromy_options = HC.MonodromyOptions(max_loops_no_progress = 10),
    start_pair_using_newton::Bool = false,
    seed = nothing,
    kwargs...,
) where {T<:Real}

    if !isnothing(seed)
        Random.seed!(seed)
    end

    variable_list = HC.variables(f)
    f_denom = generate_random_degree_2(variable_list)
    f_list = [f.expressions; f_denom]
    k = length(f_list)
    s = set_up_s(s, f)

    

    # solve critical points 
    m_1 = compute_critical_points(
        f_list,
        variable_list,
        s,
        monodromy_options,
        progress;
        start_pair_using_newton = start_pair_using_newton,
        kwargs...,
    )

    if isnothing(m_1)
        println("Couldn't solve the critical equations")
        return nothing
    end
    real_sols = real_solutions(m_1)

    logg = sum(s[i] * log(f_list[i]^2) for i = 1:k)
    ∇, Hlogg = gradient_hessian(logg, variable_list)
    ∇logg = fixed(System(∇, variables = variable_list))

    euler_char, index_list, unstable_eigenvector_lists, flag_prime =
        euler_characteristic_log(Hlogg, variable_list, real_sols)
    if flag_prime == true
        @warn "The Hessian is almost singular for some critical points"
        return nothing
    end

    set_euler_char!(progress, euler_char) # ProgressMeter

    sign_partition, failed_info = partition_of_critical_points_all_signs(
        f,
        ∇logg,
        real_sols,
        index_list,
        unstable_eigenvector_lists,
        epsilon,
        reltol,
        abstol,
        progress,
    )

    # structure output
    chambers = Vector{Chamber}()
    j = 1
    for (sign, par_list) in pairs(sign_partition)
        for par in par_list
            sub_index_list = @view index_list[par]
            μ = count_appearance(sub_index_list, length(variable_list))
            euler_char = 0
            for index in sub_index_list
                if mod(index, 2) == 1
                    euler_char -= 1
                else
                    euler_char += 1
                end
            end
            R = Chamber(sign, euler_char, μ, real_sols[par], nothing, nothing, j)
            j += 1

            push!(chambers, R)
        end
    end

    g = (f, f_denom, s)

    return ChambersResult(
        chambers,
        length(chambers),
        nnonsingular(m_1),
        length(real_sols),
        nothing,
        g,
    )
end


function count_appearance(lst, n)
    counts = zeros(Int, n + 1)
    for element in lst
        counts[element+1] += 1
    end
    return counts
end



"""
    compute_critical_points

Computes the critical points of the rational function using `HomotopyContinuation.jl`.

"""




"""
    compute_critical_points

Computes the critical points of the rational function using `HomotopyContinuation.jl`.

"""
function compute_critical_points(
    f_list::Vector{Expression},
    variable_list::Vector{Variable},
    s::Vector{T},
    monodromy_options::MonodromyOptions,
    progress::Union{Nothing,ChambersProgress};
    start_pair_using_newton::Bool = false,
    kwargs...,
) where {T<:Real}

    start_monodromy!(progress)
    show_progress = !isnothing(progress)

    k = length(f_list)
    n = length(variable_list)


    if k ≤ n
        K = n + 1
    else
        K = k
    end

    @unique_var u[1:k], v[1:(K-k)], t
    Lu = sum(u[i] * log(f_list[i]) for i = 1:k)
    Df = HC.differentiate(Lu, variable_list)
    J = HC.differentiate(Df, u)

    if start_pair_using_newton
        S = System(Df, variables = variable_list, parameters = u)
        m = HC.monodromy_solve(
            S;
            options = monodromy_options,
            show_progress = show_progress,
        )
        target_parameters = s
        if nsolutions(m) > 0
            m_1 = HC.solve(
                S,
                solutions(m);
                start_parameters = parameters(m),
                target_parameters = target_parameters,
                show_progress = show_progress,
                kwargs...,
            )

            finish_monodromy!(progress)
            return m_1
        end
    end


    A = randn(n, k)
    B = randn(n, K - k)
    Df_t = Df + t .* (A * u + B * v)
    J_t = [J + A B]

    x1 = randn(ComplexF64, n)
    p1 = compute_parameter(J_t, variable_list, x1)
    if all(abs.(p1) .> 1e-6) # first try
        S = System(Df_t, variables = variable_list, parameters = [u; v; t])
        m = HC.monodromy_solve(
            S,
            [x1],
            [p1; 1.0];
            options = monodromy_options,
            show_progress = show_progress,
            kwargs...,
        )

        target_parameters = [s; randn(ComplexF64, K - k); 0.0]
    end

    if nsolutions(m) == 0 # second try, if first failed
        S = System(Df, variables = variable_list, parameters = u)
        m = HC.monodromy_solve(
            S;
            options = monodromy_options,
            show_progress = show_progress,
            kwargs...,
        )

        finish_monodromy!(progress)
        target_parameters = s
        if nsolutions(m) == 0 # give up

            finish_monodromy!(progress)
            return nothing
        end
    end


    m_1 = HC.solve(
        S,
        solutions(m);
        start_parameters = parameters(m),
        target_parameters = target_parameters,
        show_progress = show_progress,
        kwargs...,
    )

    finish_monodromy!(progress)
    return m_1
end

function set_up_s(s, f)
    # define parameter if not provided
    if isnothing(s)
        u_last = ceil(sum(degrees(f)) / 2) + 1
        s = [ones(Float64, length(f)); -u_last]
        # check if the input parameter is correct
    else
        if length(s) != k
            error(
                "The length of the provided list of polynomials does not match the number of parameters.",
            )
        elseif any(s[1:(k-1)] .≤ 0) || s[k] ≥ 0
            error("Not all parameters have the correct sign.")
        elseif sum(degrees(f) .* s[1:(k-1)]) + 2 * s[k] > 0
            error(
                "The parameters `s` don't satisfy `-2 s_{k+1} > s_1 deg(f_1) + ... + s_k deg(f_k)``.",
            )
        end
    end

    s
end

function compute_parameter(J, variable_list, x0)
    M =
        [
            isa(entry, Expression) ? evaluate(entry, variable_list => x0) : entry for
            entry in J
        ] |> Matrix{ComplexF64}
    N = LinearAlgebra.nullspace(M)
    p = N * randn(Float64, size(N, 2))
    return p
end
