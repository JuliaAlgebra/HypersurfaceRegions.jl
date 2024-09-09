var documenterSearchIndex = {"docs":
[{"location":"#ComputingRegions.jl","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"","category":"section"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"We present a Julia package  for computing the regions (i.e., connected components) in the complement of an arrangement of real hypersurfaces in affine or projective space.","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"Our input consists of k polynomials in n variables. The output is a list of all regions of the complement of the zero set of these polynomials.","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"For each region C we report the Euler characteristic and well-chosen sample points. The list is grouped according to sign vectors sigma in  -1+1 ^k, where sigma_i is the sign of the ith polynomial on C. For hyperplane arrangements, each region is uniquely characterized by its sign vector. However, in our situation, each sign vector sigma typically corresponds to multiple connected components.","category":"page"},{"location":"#Example:-two-concentric-circles","page":"ComputingRegions.jl","title":"Example: two concentric circles","text":"","category":"section"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"Let us consider two concentric circles. For instance, we could take the two circles f_1 = x^2 + y^2 - 1=0 and f_2=x^2 + y^2 - 4=0 centered at the origin. To compute the regions of mathcalU  =    u in mathbbR^2  mid   f_1(u) cdot f_2(u)  not=  0  we can use the following code:    ","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"julia> using ComputingRegions\njulia> @var x y;\njulia> f_1 = x^2 + y^2 - 1;\njulia> f_2 = x^2 + y^2 - 4;\njulia> f = [f_1; f_2]\njulia> C = regions(f)\nRegionsResult with 3 regions:\n===============================\n9 complex critical points\n5 real critical points\n╭──────────────┬─────────────────────────────────╮\n│ sign pattern │ regions                        │\n├──────────────┼─────────────────────────────────┤\n│ + +          │ number = 1                      │\n│              │ χ = 0, μ = [1, 1, 0], unbounded │\n│ - -          │ number = 1                      │\n│              │ χ = 1, μ = [1, 0, 0], bounded   │\n│ + -          │ number = 1                      │\n│              │ χ = 0, μ = [1, 1, 0], bounded   │\n╰──────────────┴─────────────────────────────────╯","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"The output shows that mathcal U has three regions. The first region has sign pattern ++. This means that, on this region, both f_1 and f_2 are positive. On the second region, both f_1 and f_2 are negative, so it is the contractible region in the middle. The software correctly reports that this region has Euler characteristic 1. The other two regions each have one hole and thus have Euler characteristic 0. ","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"Let us visualize the three regions. We use the Plots.jl package. For this, we define a function that spits out a color corresponding to a region given a point. ","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"function region_col(x, y)\n    Ci = membership(C, [x; y]; warning = false)\n    if !isnothing(Ci)\n        return number(Ci)\n    else\n        return nothing\n    end\nend \n\nu = -4:0.05:4; v = -3:0.05:3;\nX = repeat(reshape(u, 1, :), length(v), 1);\nY = repeat(v, 1, length(u));\nZ = map(region_col, X, Y);\n\nusing Plots\ncontour(u, v, Z, fill = true, alpha = 0.25, \n                                    legend = false, \n                                    aspect_ratio = :equal)","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"This produces the following picture:","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"(Image: circles)","category":"page"},{"location":"#Example:-3-random-polynomials","page":"ComputingRegions.jl","title":"Example: 3 random polynomials","text":"","category":"section"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"We can set up a random example as follows.","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"using ComputingRegions\n@var v[1:3];\nf = [rand_poly(Float64, v, d) for d in [2, 3, 3]];\nC = regions(f)","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"Here, f consists of 3 random polynomials in v. The degrees of these polynomials are [2, 3, 3]. The coefficients are chosen from a Gaussian distribution. ","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"At a certain number of rows, the output of C in terminal is cropped. To avoid this, we can use the following command.","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"show(C; crop = false)","category":"page"},{"location":"#Documentation:-Output","page":"ComputingRegions.jl","title":"Documentation: Output","text":"","category":"section"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"RegionsResult\nRegion","category":"page"},{"location":"#ComputingRegions.RegionsResult","page":"ComputingRegions.jl","title":"ComputingRegions.RegionsResult","text":"RegionsResult\n\nA struct that collects all regions in the complement of a hypersurface arragement.\n\n\n\n\n\n","category":"type"},{"location":"#ComputingRegions.Region","page":"ComputingRegions.jl","title":"ComputingRegions.Region","text":"Region\n\nA struct that contains all information about a region.\n\n\n\n\n\n","category":"type"},{"location":"#Documentation:-Main-functions","page":"ComputingRegions.jl","title":"Documentation: Main functions","text":"","category":"section"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"regions\naffine_regions\nprojective_regions","category":"page"},{"location":"#ComputingRegions.regions","page":"ComputingRegions.jl","title":"ComputingRegions.regions","text":"regions(C::RegionsResult)\n\nReturns the vector of regions in C.\n\n\n\n\n\nregions(f::Vector{Expression})\nregions(f::System)\n\nInput a list of hypersurfaces 'f = [f1,...fk]'. Outputs the regions in the complement of the hypersurface arrangement, whether they are bounded or not, their sign patterns, Euler characteristic and the indices of the critical points in each region.\n\nOptions:\n\nδ::Float64 = 1e-4: Parameter that defines the strip around infinity.\ntarget_parameters: Specify parameters of the System f (if its has any).\nshow_progress = true: if true, prints the progress of the computation to the terminal.\nprojective_fusion = true: if true, the algorithm computes which of the regions are fused at infinity.\ns: exponents of the Morse function f_1^(s_1) * ... * f_k^(s_k) * q^(s_k+1). Here, s is a list of integers [s_1, ..., s_k, s_{k+1}] such that s_1, ..., s_k>0, s_{k+1}<0 and 2 s_{k+1} > s_1 deg(f_1) + ... + s_k deg(f_k).\nepsilon = 1e-6: how close from each critical point do we do the path tracking.\nreltol = 1e-6, abstol = 1e-9: parameters for the accuracy of the ODE solver.\nmonodromy_options = MonodromyOptions(max_loops_no_progress = 25): pass options for monodromy.\nstart_pair_using_newton::Bool = false: if true, the algorithm tries to compute a start pair for monodromy by using Newton's methods. Can reduce the number of critical points, but is less stable.\n\nExample\n\nusing ComputingRegions\n@var x y\nf = [x^2 + y^2 - 1; x^2 + y^2 - 4];\nregions(f)\n\nExample with options\n\nregions(f; δ = 1e-4, \n            monodromy_options = MonodromyOptions(max_loops_no_progress = 20))\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.affine_regions","page":"ComputingRegions.jl","title":"ComputingRegions.affine_regions","text":"affine_regions(f::Vector{Expression})\naffine_regions(f::System)\n\nInput a list of affine hypersurfaces 'f = [f1,...fk]'.  Outputs the regions in the complement of the hypersurface arrangement, their sign patterns, Euler characteristic and the indices of the critical points in each region. Accepts the same options as regions.\n\nExample\n\nusing ComputingRegions\n@var x y\nf = [x^2 + y^2 - 1; x^2 + y^2 - 4];\naffine_regions(f)\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.projective_regions","page":"ComputingRegions.jl","title":"ComputingRegions.projective_regions","text":"projective_regions(C::RegionsResult)\n\nReturns a vector of vectors. The entries of the vectors are those regions in C, which are fused in projective space.\n\nThe following code will return a vector of vectors of type Region:\n\nusing ComputingRegions\n@var x y\nf = [x^2 + y^2 - 1; x^2 + y^2 - 4];\nC = regions(f)\np = projective_regions(C)\n\n\n\n\n\n","category":"function"},{"location":"#Documentation:-Helper-functions","page":"ComputingRegions.jl","title":"Documentation: Helper functions","text":"","category":"section"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"Functions to call on Region.","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"χ\nμ\ncritical_points\nis_bounded\nis_unbounded\nis_undecided\nnumber","category":"page"},{"location":"#ComputingRegions.χ","page":"ComputingRegions.jl","title":"ComputingRegions.χ","text":"χ(C::Region)\n\nReturns the Euler characteristic.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.μ","page":"ComputingRegions.jl","title":"ComputingRegions.μ","text":"μ(C::Region)\n\nReturns the index vector.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.critical_points","page":"ComputingRegions.jl","title":"ComputingRegions.critical_points","text":"critical_points(C::Region)\n\nReturns the critical points in C.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.is_bounded","page":"ComputingRegions.jl","title":"ComputingRegions.is_bounded","text":"is_bounded(C::Region)\n\nReturns a boolean that is true, if C is bounded. \n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.is_unbounded","page":"ComputingRegions.jl","title":"ComputingRegions.is_unbounded","text":"is_unbounded(C::Region)\n\nReturns a boolean that is true, if C is undecided. \n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.is_undecided","page":"ComputingRegions.jl","title":"ComputingRegions.is_undecided","text":"is_undecided(C::Region)\n\nReturns a boolean that is true, if the algorithm could not decided whether C is bounded or not.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.number","page":"ComputingRegions.jl","title":"ComputingRegions.number","text":"number(C::Region)\n\nEach Region in a RegionsResult is assigned a number. \n\n\n\n\n\n","category":"function"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"Functions to call on RegionsResult.","category":"page"},{"location":"","page":"ComputingRegions.jl","title":"ComputingRegions.jl","text":"nregions\nnbounded\nnunbounded\nnundecided\nbounded\nunbounded\nundecided\neuler_characteristics\nindex_vectors\nncritical_complex\nncritical_real\ng\nvariables","category":"page"},{"location":"#ComputingRegions.nregions","page":"ComputingRegions.jl","title":"ComputingRegions.nregions","text":"nregions(C::RegionsResult)\n\nReturns the number of regions in F.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.nbounded","page":"ComputingRegions.jl","title":"ComputingRegions.nbounded","text":"nbounded(C::RegionsResult)\n\nReturns the number of (weakly) bounded regions in C.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.nunbounded","page":"ComputingRegions.jl","title":"ComputingRegions.nunbounded","text":"nunbounded(C::RegionsResult)\n\nReturns the number of unbounded regions in C.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.nundecided","page":"ComputingRegions.jl","title":"ComputingRegions.nundecided","text":"nundecided(C::RegionsResult)\n\nReturns the number of regions in C, where bounded or unbounded could not be decided.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.bounded","page":"ComputingRegions.jl","title":"ComputingRegions.bounded","text":"bounded(C::RegionsResult)\n\nReturns the (weakly) bounded regions in C.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.unbounded","page":"ComputingRegions.jl","title":"ComputingRegions.unbounded","text":"unbounded(C::RegionsResult)\n\nReturns the unbounded regions in C.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.undecided","page":"ComputingRegions.jl","title":"ComputingRegions.undecided","text":"undecided(C::RegionsResult)\n\nReturns the number of regions in C, where bounded or unbounded could not be decided.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.euler_characteristics","page":"ComputingRegions.jl","title":"ComputingRegions.euler_characteristics","text":"euler_characteristics(C::RegionsResult)\n\nReturns the Euler characteristics of the regions in C.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.index_vectors","page":"ComputingRegions.jl","title":"ComputingRegions.index_vectors","text":"euler_characteristics(C::RegionsResult)\n\nReturns the index vectors of the regions in C.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.ncritical_complex","page":"ComputingRegions.jl","title":"ComputingRegions.ncritical_complex","text":"ncritical_complex(C::RegionsResult)\n\nReturns the number of complex critical points of the Morse function in C.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.ncritical_real","page":"ComputingRegions.jl","title":"ComputingRegions.ncritical_real","text":"ncritical_real(C::RegionsResult)\n\nReturns the number of real critical points of the Morse function in C.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.g","page":"ComputingRegions.jl","title":"ComputingRegions.g","text":"g(C::RegionsResult)\n\nReturns the Morse function in C.\n\n\n\n\n\n","category":"function"},{"location":"#ComputingRegions.variables","page":"ComputingRegions.jl","title":"ComputingRegions.variables","text":"variables(C::RegionsResult)\n\nReturns the order of variables.\n\n\n\n\n\n","category":"function"}]
}
