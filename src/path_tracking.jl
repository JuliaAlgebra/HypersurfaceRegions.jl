
# define the gradient ascend ode
function set_up_ode(∇logg)

    function ode_log!(u, dy, y, p, t)
        n = length(∇logg)
        HC.evaluate!(u, ∇logg, @view y[:, 1])
        for i = 1:n
            dy[i] = u[i]
        end
    end

    n = length(∇logg)
    u = zeros(Float64, n)
    return (dy, y, p, t) -> ode_log!(u, dy, y, p, t)
end


# find the nearest point in a list to a given point
function finding_nearest_point(v, list_of_vectors)
    v = reshape(v, length(v))
    λ = LA.norm(v)
    distances = [LA.norm(v - x) for x in list_of_vectors]
    if minimum(distances) < 1e-4 * λ
        return argmin(distances)
    else
        return -1
    end
end



# solve ode for a given time span to find the limit critical point
function limit_critical_point_fixed_time(
    ode_log!,
    start_point::T1,
    real_sols::T2,
    time::S,
    reltol::Float64,
    abstol::Float64,
) where {T1,T2<:AbstractArray,S<:Real}

    y_0 = start_point
    tspan = (0.0, time)
    prob = DE.ODEProblem(ode_log!, y_0, tspan)
    sol = DE.solve(prob, reltol = reltol, abstol = abstol)
    y_lim = sol(time)
    critical_point_index = finding_nearest_point(y_lim, real_sols)

    return critical_point_index, sol.retcode
end


# function to find the limit critical point with various time spans 
# start from a given point
function limit_critical_point(
    ode_log!,
    start_point::T1,
    real_sols::T2,
    reltol::Float64,
    abstol::Float64,
) where {T1,T2<:AbstractArray}


    # first try t_max=10^10
    y_0 = start_point
    t_max = 10^10
    critical_point_index, retcode =
        limit_critical_point_fixed_time(ode_log!, y_0, real_sols, t_max, reltol, abstol)



    if retcode == SciMLBase.ReturnCode.Success && critical_point_index != -1
        return critical_point_index, []
    end

    # second try
    # increase time since the solution converges but does not limit to any critical point
    if retcode == SciMLBase.ReturnCode.Success && critical_point_index == -1
        t_max = 10^15
        critical_point_index, retcode =
            limit_critical_point_fixed_time(ode_log!, y_0, real_sols, t_max, reltol, abstol)
        if retcode == SciMLBase.ReturnCode.Success && critical_point_index != -1
            return critical_point_index, []
        end

        # third try
        # increase time since the solution converges but does not limit to any critical point
        if retcode == SciMLBase.ReturnCode.Success && critical_point_index == -1
            t_max = 10^20
            critical_point_index, retcode = limit_critical_point_fixed_time(
                ode_log!,
                y_0,
                real_sols,
                t_max,
                reltol,
                abstol,
            )

            if retcode == SciMLBase.ReturnCode.Success
                if critical_point_index != -1
                    return critical_point_index, []

                else
                    failed_info = [
                        y_0,
                        "does not limit to critical point for log10t_max=",
                        log10(t_max),
                    ]
                end

            else
                failed_info = ["log10t_max=", log10(t_max), "solution diverges"]
                critical_point_index = -1
            end

        else
            failed_info = ["log10t_max=", log10(t_max), "solution diverges"]
            critical_point_index = -1
        end


        # second try
        # decrease time since the solution deverges
    else
        t_max = 10^8
        critical_point_index, retcode =
            limit_critical_point_fixed_time(ode_log!, y_0, real_sols, t_max, reltol, abstol)

        if retcode == SciMLBase.ReturnCode.Success && critical_point_index != -1
            return critical_point_index, []
        end

        if retcode == SciMLBase.ReturnCode.Success && critical_point_index == -1
            failed_info = [
                "need to increase time, the solution does not limit for log10t_max=",
                log10(t_max),
            ]
        else
            # third try
            # decrease time since the solution deverges
            t_max = 10^6
            critical_point_index, retcode = limit_critical_point_fixed_time(
                ode_log!,
                y_0,
                real_sols,
                t_max,
                reltol,
                abstol,
            )
            if retcode == SciMLBase.ReturnCode.Success && critical_point_index != -1
                return critical_point_index, []
            end

            if retcode == SciMLBase.ReturnCode.Success && critical_point_index == -1
                failed_info = [
                    "need to increase time, the solution does not limit for log10t_max=",
                    log10(t_max),
                ]
            else
                failed_info = ["log10t_max=", log10(t_max), "unstable -> decrease time"]
                critical_point_index = -1
            end
        end

    end
    return critical_point_index, failed_info
end



# path tracking from a critical point to another critical point
function limit_critical_point_from_critical_point(
    ode_log!,
    real_sols::T1,
    starting_point_index::Int64,
    index_list::T2,
    unstable_eigenvector::T3,
    epsilon::Float64,
    reltol::Float64,
    abstol::Float64,
) where {T1,T2,T3<:AbstractArray}

    y_0 = real_sols[starting_point_index]
    index = index_list[starting_point_index]

    if index == 0
        return nothing, nothing
    else
        start_point = y_0 + epsilon * unstable_eigenvector
        critical_point_index, failed_info =
            limit_critical_point(ode_log!, start_point, real_sols, reltol, abstol)

        return [starting_point_index, critical_point_index], failed_info
    end
end