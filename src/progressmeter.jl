

mutable struct ChambersProgress
    progress_meter::PM.ProgressUnknown
    progress_stage::Int
    current_stage::Int
    ncritical_points::Int
    ncritical_points_classified::Int
    monodromy::Int
    euler_char::Union{Nothing,Int}
end

ChambersProgress(progress::PM.ProgressUnknown) =
    ChambersProgress(progress, 0, 0, 0, 0, 0, nothing)

# Stage 1: computing affine chambers
# Stage 2: computing projective chambers
# Stage 3: connecting critical points at infinity to affine chambers
function set_stage!(P::ChambersProgress, i::Int)
    P.current_stage = i
    update_progress!(P)
end
set_stage!(P::Nothing, i::Int) = nothing

# 0 = monodromy does not run
# 1 = monodromy runs
function start_monodromy!(P::ChambersProgress)
    P.monodromy = 1
    if P.current_stage == 1
        print(Crayon(foreground = :green), "\n\n\nComputing critical points:\n")
    elseif P.current_stage == 2
        print(
            Crayon(foreground = :green),
            "\n\n\n\nComputing critical points at infinity:\n",
        )
    end
end
start_monodromy!(P::Nothing) = nothing
function finish_monodromy!(P::ChambersProgress)
    P.monodromy = 0
    print(" \n")
end
finish_monodromy!(P::Nothing) = nothing

function set_ncritical_points!(P::ChambersProgress, i::Int)
    P.ncritical_points = i
    update_progress!(P)
end
set_ncritical_points!(P::Nothing, i::Int) = nothing

function set_ncritical_points_classified!(P::ChambersProgress, i::Int)
    P.ncritical_points_classified = i
    update_progress!(P)
end
set_ncritical_points_classified!(P::Nothing, i::Int) = nothing

function set_euler_char!(P::ChambersProgress, χ::Int)
    if P.current_stage == 1
        P.euler_char = χ
    end
    update_progress!(P)
end
set_euler_char!(P::Nothing, χ::Int) = nothing



## Show 
function showvalues(progress::ChambersProgress)

    if progress.current_stage == 0
        text = []
    elseif progress.current_stage == 1
        text = [("Computing affine chambers", "")]
        if progress.monodromy == 0 && progress.ncritical_points > 0
            n = progress.ncritical_points_classified
            N = progress.ncritical_points
            push!(text, ("Partitioning critical points into chambers", "$n/$N"))
        end
        if !isnothing(progress.euler_char)
            push!(
                text,
                ("Euler characteristic of affine chambers", "$(progress.euler_char)"),
            )
        end
    elseif progress.current_stage == 2
        text = [("Computing chambers at infinity:", "")]
        if progress.monodromy == 0 && progress.ncritical_points > 0
            n = progress.ncritical_points_classified
            N = progress.ncritical_points
            push!(text, ("Partitioning critical points into chambers", "$n/$N"))
        end
    elseif progress.current_stage == 3
        text = [("Connecting critical points at infinity to affine chambers", "")]
        if progress.ncritical_points > 0
            n = progress.ncritical_points_classified
            N = progress.ncritical_points
            push!(text, ("Critical points classified", "$n/$N"))
        end
    end


    text
end

function update_progress!(progress::ChambersProgress)
    j = progress.progress_stage + 1
    PM.update!(progress.progress_meter, j, showvalues = showvalues(progress))
end
update_progress!(progress::Nothing) = nothing

function finish_progress!(progress::ChambersProgress)
    progress.current_stage = 0
    update_progress!(progress)
    PM.finish!(progress.progress_meter, showvalues = showvalues(progress))
end
finish_progress!(progress::Nothing) = nothing
