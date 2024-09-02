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
╭──────────────┬─────────────────────────────────╮
│ sign pattern │ chambers                         │
├──────────────┼─────────────────────────────────┤
│ - -          │ number = 1                      │
│              │ χ = 1, μ = [1, 0, 0], bounded   │
│ + -          │ number = 1                      │
│              │ χ = 0, μ = [2, 2, 0], bounded   │
│ + +          │ number = 1                      │
│              │ χ = 0, μ = [2, 2, 0], unbounded │
╰──────────────┴─────────────────────────────────╯
```

The output shows that $\mathcal U$ has three chambers. The first chamber has sign pattern $++$. This means that, on this chamber, both $f_1$ and $f_2$ are positive. On the second chamber, both $f_1$ and $f_2$ are negative, so it is the contractible chamber in the middle. The software correctly reports that this chamber has Euler characteristic 1. The other two chambers each have one hole and thus have Euler characteristic 0. 

## Example: 3 random polynomials

We can set up a random example as follows.
```julia
julia> using Chambers
julia> @var v[1:3];
julia> f = [rand_poly(Float64, v, d) for d in [2, 3, 3]];
julia> C = chambers(f)
```

Here, `f` consists of 3 random polynomials in `v`. The degrees of these polynomials are `[2, 3, 3]`. The coefficients are chosen from a Gaussian distribution.

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
```