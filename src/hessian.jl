# calculate the gradient and hessian of f 
function gradient_hessian(f::Expression, variable_list::Vector{Variable})
    ∇f = HC.differentiate(f, variable_list)
    Hf = HC.differentiate(∇f, variable_list)
    return ∇f, Hf
end


# return index and unstable eigenvectors of the hessian 
function index_unstable_eigenvector(
    Hg::Array{Expression},
    variable_list::Vector{Variable},
    a::T,
) where {T<:AbstractArray}
    Hg_eval = [evaluate(h, variable_list => a) for h in Hg]
    if is_almost_singular(Hg_eval)
        flag = true
        println("The Hessian is almost singular for", a)
        return nothing, nothing, flag
    else
        flag = false
    end
    eigvals, eigvecs = LA.eigen!(Hg_eval)
    index = sum(eigvals .> 0)
    unstable_eigenvectors = eigvecs[:, eigvals.>0]
    return index, unstable_eigenvectors, flag
end


## calculate the index, unstable eigenvectors of the critical points and the Euler characteristic of the hypersurface arrangement
function euler_characteristic_log(
    Hf::Array{Expression},
    variable_list::Vector{Variable},
    real_sols::T,
) where {T<:AbstractArray}
    index_list::Vector{Int} = []
    unstable_vector_list::Vector{Matrix{Float64}} = []
    flag_prime = false
    for i = 1:length(real_sols)
        index, unstable_eigenvectors, flag =
            index_unstable_eigenvector(Hf, variable_list, real_sols[i])
        if flag == true
            flag_prime = true
            return nothing, nothing, nothing, flag_prime
        end

        push!(index_list, index)
        push!(unstable_vector_list, unstable_eigenvectors)
    end

    euler_char = count(index_list .== 0) - count(index_list .== 1)
    return euler_char, index_list, unstable_vector_list, flag_prime
end


# test if a matrix is almost singular
function is_almost_singular(matrix::Matrix{Float64}, threshold = 1e10)
    condition_number = LA.cond(matrix)
    return condition_number > threshold
end
