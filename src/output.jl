export RegionsResult,
    Region,
    χ,
    μ,
    critical_points,
    is_bounded,
    g,
    regions,
    nregions,
    nbounded,
    nunbounded,
    nundecided,
    bounded,
    unbounded,
    undecided,
    is_bounded,
    is_unbounded,
    is_undecided,
    euler_characteristics,
    index_vectors,
    ncritical_complex,
    ncritical_real,
    projective_regions,
    sign,
    number,
    variables

"""
    Region

A struct that contains all information about a region.
"""
struct Region
    sign::Vector{Int}
    χ::Int
    μ::Vector{Int}
    critical_points::Union{Vector{Float64},Vector{Vector{Float64}}}
    g::Union{Nothing,Tuple{System,Expression,Vector{Int}}}
    is_bounded::Union{Nothing,Int} # 0 = unbounded, 1 = bounded, 2 = undecided
    region_number::Int
end

"""
    sign(C::Region)

Returns the sign vector.
"""
Base.sign(C::Region) = C.sign

"""
    χ(C::Region)

Returns the Euler characteristic.
"""
χ(C::Region) = C.χ

"""
    μ(C::Region)

Returns the index vector.
"""
μ(C::Region) = C.μ

"""
    critical_points(C::Region)

Returns the critical points in `C`.
"""
critical_points(C::Region) = C.critical_points

"""
    is_unbounded(C::Region)

Returns a boolean that is `true`, if `C` is undecided. 
"""
is_unbounded(C::Region) = C.is_bounded == 0

"""
    is_bounded(C::Region)

Returns a boolean that is `true`, if `C` is bounded. 
"""
is_bounded(C::Region) = C.is_bounded == 1

"""
    is_undecided(C::Region)

Returns a boolean that is `true`, if the algorithm could not decided whether `C` is bounded or not.
"""
is_undecided(C::Region) = C.is_bounded == 2

"""
    number(C::Region)

Each `Region` in a `RegionsResult` is assigned a number. 
"""
number(C::Region) = C.region_number

g(C::Region) = C.g


"""
    RegionsResult

A struct that collects all regions in the complement of a hypersurface arragement.
"""
struct RegionsResult
    region_list::Vector{Region}
    nregions::Int
    ncritical_complex::Int
    ncritical_real::Int
    projective_regions::Union{Nothing,Vector{Vector{Int}}}
    g::Tuple{System,Expression,Vector{Int}}
end

"""
    regions(C::RegionsResult)

Returns the vector of regions in `C`.
"""
function regions(C::RegionsResult)
    out = map(C.region_list) do Cᵢ
        Region(Cᵢ.sign, Cᵢ.χ, Cᵢ.μ, Cᵢ.critical_points, C.g, Cᵢ.is_bounded, Cᵢ.region_number)
    end
    return out
end

"""
    projective_regions(C::RegionsResult)

Returns a vector of vectors. The entries of the vectors are those regions in `C`, which are fused in projective space.

The following code will return a vector of vectors of type `Region`:
```julia
using HypersurfaceRegions
@var x y
f = [x^2 + y^2 - 1; x^2 + y^2 - 4];
C = regions(f)
p = projective_regions(C)
```
"""
function projective_regions(C::RegionsResult)
    c = regions(C)
    p = C.projective_regions
    P = map(p) do pᵢ
        filter(ci -> in(number(ci), pᵢ), c)
    end

    P
end

"""
    nregions(C::RegionsResult)

Returns the number of regions in F.
"""
nregions(C::RegionsResult) = C.nregions

"""
    nbounded(C::RegionsResult)

Returns the number of (weakly) bounded regions in `C`.
"""
nbounded(C::RegionsResult) = count(is_bounded, regions(C))

"""
    nunbounded(C::RegionsResult)

Returns the number of unbounded regions in `C`.
"""
nunbounded(C::RegionsResult) = count(is_unbounded, regions(C))

"""
    nundecided(C::RegionsResult)

Returns the number of regions in `C`, where bounded or unbounded could not be decided.
"""
nundecided(C::RegionsResult) = count(is_undecided, regions(C))

"""
    bounded(C::RegionsResult)

Returns the (weakly) bounded regions in `C`.
"""
bounded(C::RegionsResult) = filter(is_bounded, regions(C))

"""
    unbounded(C::RegionsResult)

Returns the unbounded regions in `C`.
"""
unbounded(C::RegionsResult) = filter(is_unbounded, regions(C))

"""
    undecided(C::RegionsResult)

Returns the number of regions in `C`, where bounded or unbounded could not be decided.
"""
undecided(C::RegionsResult) = filter(is_undecided, regions(C))

"""
    euler_characteristics(C::RegionsResult)

Returns the Euler characteristics of the regions in `C`.
"""
euler_characteristics(C::RegionsResult) = map(r -> χ(r), regions(C))
χ(C::RegionsResult) = euler_characteristics(C)

"""
    euler_characteristics(C::RegionsResult)

Returns the index vectors of the regions in `C`.
"""
index_vectors(C::RegionsResult) = map(r -> μ(r), regions(C))
μ(C::RegionsResult) = index_vectors(C)

"""
    ncritical_complex(C::RegionsResult)

Returns the number of complex critical points of the Morse function in `C`.
"""
ncritical_complex(C::RegionsResult) = C.ncritical_complex

"""
    ncritical_real(C::RegionsResult)

Returns the number of real critical points of the Morse function in `C`.
"""
ncritical_real(C::RegionsResult) = C.ncritical_real

"""
    g(C::RegionsResult)

Returns the Morse function in `C`.
"""
g(C::RegionsResult) = C.g


"""
    variables(C::RegionsResult)

Returns the order of variables.
"""
function variables(C::RegionsResult)
    f, _, _ = g(C)
    HC.HC.variables(f)
end

###############
### Show ###
###############
Base.show(C::RegionsResult; crop = true) = Base.show(stdout, C, crop = crop)
function Base.show(io::IO, C::RegionsResult; crop = true)

    header = "RegionsResult with $(C.nregions) regions:"
    println(io, header)
    println(io, "="^(length(header)))

    println(io, "$(C.ncritical_complex) complex critical points")
    println(io, "$(C.ncritical_real) real critical points")

    all_regions = C.region_list
    sign_list = [region.sign for region in all_regions]
    unique_signs = unique(sign_list)

    k = length(unique_signs) + nregions(C)
    table = Matrix{String}(undef, k, 2)

    which_sign = [findall(sign -> sign == s, sign_list) for s in unique_signs]

    i = 1
    for (w, s) in zip(which_sign, unique_signs)
        table[i, 1] = join([v == 1 ? "+" : v == -1 ? "-" : string(v) for v in s], " ")
        table[i, 2] = "number = $(length(w))"
        i += 1

        for wᵢ in w
            region = all_regions[wᵢ]
            sign = region.sign
            b = region.is_bounded
            χ = region.χ
            μ = region.μ

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
        crayon = crayon"magenta",
    )

    pretty_table(
        table;
        header = ["sign pattern", "regions"],
        header_crayon = crayon"green",
        tf = tf_unicode_rounded,
        alignment = :l,
        display_size = ds,
        highlighters = (h1, h2, h3),
    )
end

function Base.show(io::IO, C::Region)
    sign = C.sign
    sign_string = join([v == 1 ? "+" : v == -1 ? "-" : string(v) for v in sign], " ")

    header = "Region with sign pattern ($sign_string) :"
    println(io, header)

    b = C.is_bounded
    χ = C.χ
    μ = C.μ

    if isnothing(b)
        println(io, " χ = $χ, μ = $μ")
    elseif b == 0
        println(io, "χ = $χ, μ = $μ, unbounded")
    elseif b == 1
        println(io, "χ = $χ, μ = $μ, bounded")
    elseif b == 2
        println(io, "χ = $χ, μ = $μ, undecided")
    end
end
