# Chambers.jl

We present a Julia package 
for computing the connected components in the complement of an arrangement of real
hypersurfaces in affine or projective space.

Our input consists of
$k$ polynomials in $n$ variables.
The output is a list of all *chambers* (i.e. connected components) of the complement of the zero set of these polynomials.

For each chamber $C$ we report the Euler characteristic and well-chosen
sample points.
The list is grouped according to sign vectors $\sigma \in  \{-1,+1 \}^k$,
where $\sigma_i$ is the sign of the $i$th polynomial on $C$.
For hyperplane arrangements, each chamber is
uniquely characterized by its sign vector. However, in our situation, each sign vector $\sigma$ typically
corresponds to multiple connected components.

## Example: two concentric circles

Let us consider two concentric circles. For instance, we could take the two circles $f_1 = x^2 + y^2 - 1=0$ and $f_2=x^2 + y^2 - 4=0$ centered at the origin. To compute the chambers of $\mathcal{U}  =   \{ u \in \mathbb{R}^2  \mid   f_1(u) \cdot f_2(u)  \not=  0 \}$ we can use the following code:    

```julia
julia> using Chambers
julia> @var x y;
julia> f_1 = x^2 + y^2 - 1;
julia> f_2 = x^2 + y^2 - 4;
julia> f = [f_1; f_2]
julia> C = chambers(f)
ChambersResult with 3 chambers:
=============================
9 complex critical points
9 real critical points
╭──────────────┬──────────────────────────────────────────╮
│ sign pattern │ chambers                                 │
├──────────────┼──────────────────────────────────────────┤
│ - -          │ number = 1                               │
│              │ χ = 1, μ = [1, 0, 0], (weakly) bounded   │
│ + -          │ number = 1                               │
│              │ χ = 0, μ = [2, 2, 0], (weakly) bounded   │
│ + +          │ number = 1                               │
│              │ χ = 0, μ = [2, 2, 0], unbounded          │
╰──────────────┴──────────────────────────────────────────╯
```

The output shows that $\mathcal U$ has three chambers. The first chamber has sign pattern $++$. This means that, on this chamber, both $f_1$ and $f_2$ are positive. On the second chamber, both $f_1$ and $f_2$ are negative, so it is the contractible chamber in the middle. The software correctly reports that this chamber has Euler characteristic 1. The other two chambers each have one hole and thus have Euler characteristic 0. 

Let us visualize the three chambers. We use the [Plots.jl](https://docs.juliaplots.org/) package. For this, we define a function that spits out a color corresponding to a region given a point. It uses the [membership](@ref) function. 

```julia
function chamber_col(x, y)
    Ci = membership(C, [x; y]; warning = false)
    if !isnothing(Ci)
        return number(Ci)
    else
        return nothing
    end
end 

u = -4:0.05:4; v = -3:0.05:3;
X = repeat(reshape(u, 1, :), length(v), 1);
Y = repeat(v, 1, length(u));
Z = map(chamber_col, X, Y);

using Plots
contour(u, v, Z, fill = true, alpha = 0.25, 
                                    legend = false, 
                                    aspect_ratio = :equal)
```

This produces the following picture:

![circles](circles.png)

## Example: 3 random polynomials

We can set up a random example as follows.
```julia
using Chambers
@var v[1:3];
f = [rand_poly(Float64, v, d) for d in [2, 3, 3]];
C = chambers(f)
```

Here, `f` consists of 3 random polynomials in `v`. The degrees of these polynomials are `[2, 3, 3]`. The coefficients are chosen from a Gaussian distribution. 

At a certain number of rows, the output of `C` in terminal is cropped. To avoid this, we can use the following command.
```julia
show(C; crop = false)
```


## Documentation: Output

```@docs
ChambersResult
Chamber
```

## Documentation: Main functions

```@docs
chambers
affine_chambers
projective_chambers
```


## Documentation: Helper functions

Functions to call on [Chamber](@ref).
```@docs
χ
μ
critical_points
is_bounded
number
```

Functions to call on [ChambersResult](@ref).
```@docs
nchambers
nbounded
nunbounded
bounded
unbounded
euler_characteristics
index_vectors
ncritical_complex
ncritical_real
g
variables
```