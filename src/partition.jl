


# input a list, output a list of indices that have the same value in the list
function partition_indices(lst)
    partitions = Dict()
    for (index, element) in enumerate(lst)
        if haskey(partitions, element)
            push!(partitions[element], index)
        else
            partitions[element] = [index]
        end
    end
    return partitions
end


function partition_points_via_sign(f::System, real_sols::Vector{Vector{Float64}})
    sign_list_real_sols = map(real_sols) do real_sol
        Base.sign.(f(real_sol))
    end
    return partition_indices(sign_list_real_sols)
end




function generate_random_degree_2(variable_list::Vector{Variable})
    c_1 = randn(length(variable_list))
    c_2 = randn(length(variable_list))
    c_3 = 1
    l = 0
    g = 0
    for (i, variable) in enumerate(variable_list)
        l += c_1[i] * variable
        g += (c_2[i] - variable)^2
    end
    g += l^2 + c_3^2
    return g
end




function partition_of_critical_points_single_sign(
    ode_log!,
    real_sols::Vector{Vector{Float64}},
    critical_points_indices::Vector{Int},
    index_list::Vector{Int},
    unstable_eigenvector_list,
    epsilon::Float64,
    reltol::Float64,
    abstol::Float64,
)


    sub_index_list = @view index_list[critical_points_indices]
    sub_real_sols = @view real_sols[critical_points_indices]
    graph = LG.SimpleGraph(length(critical_points_indices))
    connectivity_status = zeros(Int, length(critical_points_indices))

    # count the number of index 0 critical points
    count_index_0 = sum([index == 0 for index in sub_index_list])



    if count_index_0 == 1
        # we do not need to do any path tracking in this case
        index_0_index = in(0, sub_index_list)
        return [critical_points_indices], []
    end


    failed_info_list = []
    for i = 1:length(sub_index_list)
        if connectivity_status[i] == 0 && sub_index_list[i] == 1
            # need to do path tracking in two directions
            critical_point_index = critical_points_indices[i]
            unstable_eigenvector = unstable_eigenvector_list[critical_point_index]
            pair_pos, failed_info_pos = limit_critical_point_from_critical_point(
                ode_log!,
                sub_real_sols,
                i,
                sub_index_list,
                unstable_eigenvector,
                epsilon,
                reltol,
                abstol,
            )
            if isempty(failed_info_pos)
                LG.add_edge!(graph, pair_pos[1], pair_pos[2])
            else
                push!(
                    failed_info_list,
                    [critical_point_index, epsilon, unstable_eigenvector, failed_info_pos],
                )
            end

            pair_neg, failed_info_neg = limit_critical_point_from_critical_point(
                ode_log!,
                sub_real_sols,
                i,
                sub_index_list,
                unstable_eigenvector,
                -epsilon,
                reltol,
                abstol,
            )

            if isempty(failed_info_neg)
                LG.add_edge!(graph, pair_neg[1], pair_neg[2])
            else
                push!(
                    failed_info_list,
                    [critical_point_index, -epsilon, unstable_eigenvector, failed_info_neg],
                )
            end
            connectivity_status[i] = 1
        end
    end


    for i = 1:length(sub_index_list)
        if connectivity_status[i] == 0 && sub_index_list[i] > 1
            # track paths and stop whenever one path converges
            critical_point_index = critical_points_indices[i]
            sub_failed_info_list = []
            for v in eachcol(unstable_eigenvector_list[critical_point_index])
                pair, failed_info = limit_critical_point_from_critical_point(
                    ode_log!,
                    sub_real_sols,
                    i,
                    sub_index_list,
                    v,
                    epsilon,
                    reltol,
                    abstol,
                )
                if isempty(failed_info)
                    LG.add_edge!(graph, pair[1], pair[2])
                    connectivity_status[i] = 1
                    sub_failed_info_list = []
                    break
                else
                    push!(
                        sub_failed_info_list,
                        [critical_point_index, epsilon, v, failed_info],
                    )
                end

                pair, failed_info = limit_critical_point_from_critical_point(
                    ode_log!,
                    sub_real_sols,
                    i,
                    sub_index_list,
                    v,
                    -epsilon,
                    reltol,
                    abstol,
                )
                if isempty(failed_info)
                    LG.add_edge!(graph, pair[1], pair[2])
                    connectivity_status[i] = 1
                    sub_failed_info_list = []
                    break
                else
                    push!(
                        sub_failed_info_list,
                        [critical_point_index, epsilon, v, failed_info],
                    )
                end
            end

            if !isempty(sub_failed_info_list)
                push!(failed_info_list, sub_failed_info_list)
            end
        end
    end


    partition = LG.connected_components(graph)
    partition_critical_point_indices = []
    for par in partition
        push!(partition_critical_point_indices, @view(critical_points_indices[par]))
    end
    return partition_critical_point_indices, failed_info_list
end



function partition_of_critical_points_all_signs(
    f::System,
    ∇logg::T,
    real_sols::Vector{Vector{Float64}},
    index_list::Vector{Int},
    unstable_eigenvector_list,
    epsilon::Float64,
    reltol::Float64,
    abstol::Float64,
    progress::Union{Nothing,ChambersProgress},
) where {T<:AbstractSystem}

    ode_log! = set_up_ode(∇logg)

    sign_critical_point_indices_dict = partition_points_via_sign(f, real_sols)
    sign_partition = Dict()
    failed_info = []

    P = pairs(sign_critical_point_indices_dict)

    N = length(real_sols)
    i = 0
    set_ncritical_points!(progress, N)
    set_ncritical_points_classified!(progress, i)

    for (sign_pattern, critical_point_indices) in P
        single_sign_par, failed_info_list = partition_of_critical_points_single_sign(
            ode_log!,
            real_sols,
            critical_point_indices,
            index_list,
            unstable_eigenvector_list,
            epsilon,
            reltol,
            abstol,
        )
        sign_partition[sign_pattern] = single_sign_par
        push!(failed_info, failed_info_list)

        i += length(critical_point_indices)
        set_ncritical_points_classified!(progress, i)
    end

    return sign_partition, failed_info
end
