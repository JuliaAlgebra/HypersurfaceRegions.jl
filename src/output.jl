export ChambersResult,
    Chamber,
    χ,
    μ,
    critical_points,
    is_bounded,
    g,
    chambers,
    nchambers,
    nbounded,
    nunbounded,
    bounded,
    unbounded,
    euler_characteristics,
    index_vectors,
    ncritical_complex,
    ncritical_real,
    projective_chambers,
    sign,
    number,
    variables

"""
    Chamber

A struct that contains all information about a chamber.
"""
struct Chamber
    sign::Vector{Int}
    χ::Int
    μ::Vector{Int}
    critical_points::Union{Vector{Float64},Vector{Vector{Float64}}}
    g::Union{Nothing,Tuple{System,Expression,Vector{Int}}}
    is_bounded::Union{Nothing,Int} # 0 = unbounded, 1 = bounded, 2 = undecided
    chamber_number::Int
end

"""
    sign(C::Chamber)

Returns the sign vector.
"""
Base.sign(C::Chamber) = C.sign

"""
    χ(C::Chamber)

Returns the Euler characteristic.
"""
χ(C::Chamber) = C.χ

"""
    μ(C::Chamber)

Returns the index vector.
"""
μ(C::Chamber) = C.μ

"""
    critical_points(C::Chamber)

Returns the critical points in `C`.
"""
critical_points(C::Chamber) = C.critical_points

"""
    is_unbounded(C::Chamber)

Returns a boolean that is `true`, if `C` is undecided. 
"""
is_bounded(C::Chamber) = C.is_bounded == 0

"""
    is_bounded(C::Chamber)

Returns a boolean that is `true`, if `C` is bounded. 
"""
is_bounded(C::Chamber) = C.is_bounded == 1

"""
    is_undecided(C::Chamber)

Returns a boolean that is `true`, if the algorithm could not decided whether `C` is bounded or not.
"""
is_bounded(C::Chamber) = C.is_bounded == 2

"""
    number(C::Chamber)

Each `Chamber` in a `ChambersResult` is assigned a number. 
"""
number(C::Chamber) = C.chamber_number

g(C::Chamber) = C.g


"""
    ChambersResult

A struct that collects all chambers in the complement of a hypersurface arragement.
"""
struct ChambersResult
    chamber_list::Vector{Chamber}
    nchambers::Int
    ncritical_complex::Int
    ncritical_real::Int
    projective_chambers::Union{Nothing,Vector{Vector{Int}}}
    g::Tuple{System,Expression,Vector{Int}}
end

"""
    chambers(C::ChambersResult)

Returns the vector of chambers in `C`.
"""
function chambers(C::ChambersResult)
    out = map(C.chamber_list) do Cᵢ
        Chamber(Cᵢ.sign, Cᵢ.χ, Cᵢ.μ, Cᵢ.critical_points, C.g, Cᵢ.is_bounded, Cᵢ.chamber_number)
    end
    return out
end

"""
    projective_chambers(C::ChambersResult)

Returns a vector of vectors. The entries of the vectors are those chambers in `C`, which are fused in projective space.

The following code will return a vector of vectors of type `Chamber`:
```julia
using Chambers
@var x y
f = [x^2 + y^2 - 1; x^2 + y^2 - 4];
C = chambers(f)
p = projective_chambers(C)
```
"""
function projective_chambers(C::ChambersResult) 
    c = chambers(C)
    p = C.projective_chambers
    P = map(p) do pᵢ
        filter(ci -> in(number(ci), pᵢ), c)
    end

    P
end

"""
    nchambers(C::ChambersResult)

Returns the number of chambers in F.
"""
nchambers(C::ChambersResult) = C.nchambers

"""
    nbounded(C::ChambersResult)

Returns the number of (weakly) bounded chambers in `C`.
"""
nbounded(C::ChambersResult) = count(is_bounded, chambers(C))

"""
    nunbounded(C::ChambersResult)

Returns the number of unbounded chambers in `C`.
"""
nunbounded(C::ChambersResult) = count(is_unbounded, chambers(C))

"""
    nundecided(C::ChambersResult)

Returns the number of chambers in `C`, where bounded or unbounded could not be decided.
"""
nunbounded(C::ChambersResult) = count(is_undecided, chambers(C))

"""
    bounded(C::ChambersResult)

Returns the (weakly) bounded chambers in `C`.
"""
bounded(C::ChambersResult) = filter(is_bounded, chambers(C))

"""
    unbounded(C::ChambersResult)

Returns the unbounded chambers in `C`.
"""
unbounded(C::ChambersResult) = filter(is_unbounded, chambers(C))

"""
    undecided(C::ChambersResult)

Returns the number of chambers in `C`, where bounded or unbounded could not be decided.
"""
nunbounded(C::ChambersResult) = filter(is_undecided, chambers(C))

"""
    euler_characteristics(C::ChambersResult)

Returns the Euler characteristics of the chambers in `C`.
"""
euler_characteristics(C::ChambersResult) = map(r -> χ(r), chambers(C))
χ(C::ChambersResult) = euler_characteristics(C)

"""
    euler_characteristics(C::ChambersResult)

Returns the index vectors of the chambers in `C`.
"""
index_vectors(C::ChambersResult) = map(r -> μ(r), chambers(C))
μ(C::ChambersResult) = index_vectors(C)

"""
    ncritical_complex(C::ChambersResult)

Returns the number of complex critical points of the Morse function in `C`.
"""
ncritical_complex(C::ChambersResult) = C.ncritical_complex

"""
    ncritical_real(C::ChambersResult)

Returns the number of real critical points of the Morse function in `C`.
"""
ncritical_real(C::ChambersResult) = C.ncritical_real

"""
    g(C::ChambersResult)

Returns the Morse function in `C`.
"""
g(C::ChambersResult) = C.g


"""
    variables(C::ChambersResult)

Returns the order of variables.
"""
function variables(C::ChambersResult)
    f, _, _ = g(C)
    HC.HC.variables(f)
end

###############
### Show ###
###############
Base.show(C::ChambersResult; crop = true) = Base.show(stdout, R, crop = crop)
function Base.show(io::IO, C::ChambersResult; crop = true)

    header = "ChambersResult with $(C.nchambers) chambers:"
    println(io, header)
    println(io, "="^(length(header)))

    println(io, "$(C.ncritical_complex) complex critical points")
    println(io, "$(C.ncritical_real) real critical points")

    all_chambers = C.chamber_list
    sign_list = [chamber.sign for chamber in all_chambers]
    unique_signs = unique(sign_list)

    k = length(unique_signs) + nchambers(C)
    table = Matrix{String}(undef, k, 2)

    which_sign = [findall(sign -> sign == s, sign_list) for s in unique_signs]

    i = 1
    for (w, s) in zip(which_sign, unique_signs)
        table[i, 1] = join([v == 1 ? "+" : v == -1 ? "-" : string(v) for v in s], " ")
        table[i, 2] = "number = $(length(w))"
        i += 1
        for wᵢ in w
            chamber = all_chambers[wᵢ]
            sign = chamber.sign
            b = chamber.is_bounded
            χ = chamber.χ
            μ = chamber.μ

            table[i, 1] = ""
            if isnothing(b)
                table[i, 2] = " χ = $χ, μ = $μ"
            elseif b == 0
                table[i, 2] = "χ = $χ, μ = $μ, unbounded"
            elseif b == 1
                table[i, 2] = "χ = $χ, μ = $μ, bounded"
            elseif b == 2
                table[i, 2] = "χ = $χ, μ = $μ, undecided"
            end
            i += 1
        end
    end

    ds = displaysize()
    if !crop
        ds = (k + 8, last(ds))
    end


    h1 = Highlighter(
        f = (data, i, j) -> (last(data[i, j], 8) == " bounded");
        crayon = crayon"208",
    )
    h2 = Highlighter(
        f = (data, i, j) -> (last(data[i, j], 9) == "unbounded");
        crayon = crayon"blue",
    )
    h3 = Highlighter(
        f = (data, i, j) -> (last(data[i, j], 9) == "undecided");
        crayon = crayon"cyan",
    )

    pretty_table(
        table;
        header = ["sign pattern", "chambers"],
        header_crayon = crayon"green",
        tf = tf_unicode_rounded,
        alignment = :l,
        display_size = ds,
        highlighters = (h1, h2),
    )
end

function Base.show(io::IO, C::Chamber)
    sign = C.sign
    sign_string = join([v == 1 ? "+" : v == -1 ? "-" : string(v) for v in sign], " ")

    header = "Chamber with sign pattern ($sign_string) :"
    println(io, header)

    b = C.is_bounded
    χ = C.χ
    μ = C.μ

    if isnothing(b)
        table[i, 2] = " χ = $χ, μ = $μ"
    elseif b == 0
        table[i, 2] = "χ = $χ, μ = $μ, unbounded"
    elseif b == 1
        table[i, 2] = "χ = $χ, μ = $μ, bounded"
    elseif b == 2
        table[i, 2] = "χ = $χ, μ = $μ, undecided"
    end
end
