
[![][docs-stable-img]][docs-stable-url]

# HypersurfaceRegions.jl

We present a Julia package 
for computing the *regions* (i.e., connected components) in the complement of an arrangement of real
hypersurfaces.

Our input consists of
$k$ polynomials in $n$ variables:

$$
f_1,f_2,\ldots, f_k  \in  \mathbb{R}[x_1,\ldots,x_n]. 
$$

The output is a list of all regions of the $n$-dimensional manifold

$$\mathcal{U}    =     \\{   u \in \mathbb{R}^n  \mid  f_1(u) \cdot f_2(u)  \cdots  f_k(u)  \not=  0   \\}. 
$$

For each region $C$ we report the Euler characteristic and well-chosen
sample points.
The list is grouped according to sign vectors $\sigma \in  \\{-1,+1 \\}^k$,
where $\sigma_i$ is the sign of $f_i$ on $C$.
For hyperplane arrangements, each region is
uniquely characterized by its sign vector. However, in our situation, each sign vector $\sigma$ typically
corresponds to multiple connected components.

## Example 

Let us consider two concentric circles. For instance, we could take the two circles $f_1 = x^2 + y^2 - 1=0$ and $f_2=x^2 + y^2 - 4=0$ centered at the origin. To compute the regions of $\mathcal{U}  =   \\{ u \in \mathbb{R}^2  \mid   f_1(u) \cdot f_2(u)  \not=  0 \\}$ we can use the following code:    

```julia
julia> using HypersurfaceRegions
julia> @var x y;
julia> f_1 = x^2 + y^2 - 1;
julia> f_2 = x^2 + y^2 - 4;
julia> f = [f_1; f_2]
julia> regions(f)
RegionsResult with 3 regions:
===============================
9 complex critical points
5 real critical points
╭──────────────┬─────────────────────────────────╮
│ sign pattern │ regions                        │
├──────────────┼─────────────────────────────────┤
│ + +          │ number = 1                      │
│              │ χ = 0, μ = [1, 1, 0], unbounded │
│ - -          │ number = 1                      │
│              │ χ = 1, μ = [1, 0, 0], bounded   │
│ + -          │ number = 1                      │
│              │ χ = 0, μ = [1, 1, 0], bounded   │
╰──────────────┴─────────────────────────────────╯
```

The output shows that $\mathcal U$ has three regions. The first region has sign pattern $++$. This means that, on this region, both $f_1$ and $f_2$ are positive. On the second region, both $f_1$ and $f_2$ are negative, so it is the contractible region in the middle. The software correctly reports that this region has Euler characteristic 1. The other two regions each have one hole and thus have Euler characteristic 0. 

[docs-stable-img]: https://img.shields.io/badge/docs-online-blue.svg
[docs-stable-url]: https://juliaalgebra.github.io/HypersurfaceRegions.jl/
