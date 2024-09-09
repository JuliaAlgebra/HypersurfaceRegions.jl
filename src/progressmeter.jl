

mutable struct RegionsProgress
    progress_meter::PM.ProgressUnknown
    progress_stage::Int
    current_stage::Int
    ncritical_points::Int
    ncritical_points_classified::Int
    monodromy::Int
    euler_char::Union{Nothing,Int}
end

RegionsProgress(progress::PM.ProgressUnknown) =
    RegionsProgress(progress, 0, 0, 0, 0, 0, nothing)

# Stage 1: computing affine regions
# Stage 2: computing projective regions
# Stage 3: connecting critical points at infinity to affine regions
function set_stage!(P::RegionsProgress, i::Int)
    P.current_stage = i
    update_progress!(P)
end
set_stage!(P::Nothing, i::Int) = nothing

# 0 = monodromy does not run
# 1 = monodromy runs
function start_monodromy!(P::RegionsProgress)
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
function finish_monodromy!(P::RegionsProgress)
    P.monodromy = 0
    print(" \n")
end
finish_monodromy!(P::Nothing) = nothing

function set_ncritical_points!(P::RegionsProgress, i::Int)
    P.ncritical_points = i
    update_progress!(P)
end
set_ncritical_points!(P::Nothing, i::Int) = nothing

function set_ncritical_points_classified!(P::RegionsProgress, i::Int)
    P.ncritical_points_classified = i
    update_progress!(P)
end
set_ncritical_points_classified!(P::Nothing, i::Int) = nothing

function set_euler_char!(P::RegionsProgress, χ::Int)
    if P.current_stage == 1
        P.euler_char = χ
    end
    update_progress!(P)
end
set_euler_char!(P::Nothing, χ::Int) = nothing



## Show 
function showvalues(progress::RegionsProgress)

    if progress.current_stage == 0
        text = []
    elseif progress.current_stage == 1
        text = [("Computing affine regions", "")]
        if progress.monodromy == 0 && progress.ncritical_points > 0
            n = progress.ncritical_points_classified
            N = progress.ncritical_points
            push!(text, ("Partitioning critical points into regions", "$n/$N"))
        end
        if !isnothing(progress.euler_char)
            push!(
                text,
                ("Euler characteristic of affine regions", "$(progress.euler_char)"),
            )
        end
    elseif progress.current_stage == 2
        text = [("Computing regions at infinity:", "")]
        if progress.monodromy == 0 && progress.ncritical_points > 0
            n = progress.ncritical_points_classified
            N = progress.ncritical_points
            push!(text, ("Partitioning critical points into regions", "$n/$N"))
        end
    elseif progress.current_stage == 3
        text = [("Connecting critical points at infinity to affine regions", "")]
        if progress.ncritical_points > 0
            n = progress.ncritical_points_classified
            N = progress.ncritical_points
            push!(text, ("Critical points classified", "$n/$N"))
        end
    end


    text
end

function update_progress!(progress::RegionsProgress)
    j = progress.progress_stage + 1
    PM.update!(progress.progress_meter, j, showvalues = showvalues(progress))
end
update_progress!(progress::Nothing) = nothing

function finish_progress!(progress::RegionsProgress)
    progress.current_stage = 0
    update_progress!(progress)
    PM.finish!(progress.progress_meter, showvalues = showvalues(progress))
end
finish_progress!(progress::Nothing) = nothing
