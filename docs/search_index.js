var documenterSearchIndex = {"docs":
[{"location":"#Chambers.jl","page":"Chambers.jl","title":"Chambers.jl","text":"","category":"section"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"We present a Julia package  for computing the connected components in the complement of an arrangement of real hypersurfaces in affine or projective space.","category":"page"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"Our input consists of k polynomials in n variables. The output is a list of all chambers (i.e. connected components) of the complement of the zero set of these polynomials.","category":"page"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"For each chamber C we report the Euler characteristic and well-chosen sample points. The list is grouped according to sign vectors sigma in  -1+1 ^k, where sigma_i is the sign of the ith polynomial on C. For hyperplane arrangements, each chamber is uniquely characterized by its sign vector. However, in our situation, each sign vector sigma typically corresponds to multiple connected components.","category":"page"},{"location":"#Example:-two-concentric-circles","page":"Chambers.jl","title":"Example: two concentric circles","text":"","category":"section"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"Let us consider two concentric circles. For instance, we could take the two circles f_1 = x^2 + y^2 - 1=0 and f_2=x^2 + y^2 - 4=0 centered at the origin. To compute the chambers of mathcalU  =    u in mathbbR^2  mid   f_1(u) cdot f_2(u)  not=  0  we can use the following code:    ","category":"page"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"julia> using Chambers\njulia> @var x y;\njulia> f_1 = x^2 + y^2 - 1;\njulia> f_2 = x^2 + y^2 - 4;\njulia> f = [f_1; f_2]\njulia> C = chambers(f)\nChambersResult with 3 chambers:\n=============================\n9 complex critical points\n9 real critical points\n╭──────────────┬─────────────────────────────────╮\n│ sign pattern │ chambers                         │\n├──────────────┼─────────────────────────────────┤\n│ - -          │ number = 1                      │\n│              │ χ = 1, μ = [1, 0, 0], bounded   │\n│ + -          │ number = 1                      │\n│              │ χ = 0, μ = [2, 2, 0], bounded   │\n│ + +          │ number = 1                      │\n│              │ χ = 0, μ = [2, 2, 0], unbounded │\n╰──────────────┴─────────────────────────────────╯","category":"page"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"The output shows that mathcal U has three chambers. The first chamber has sign pattern ++. This means that, on this chamber, both f_1 and f_2 are positive. On the second chamber, both f_1 and f_2 are negative, so it is the contractible chamber in the middle. The software correctly reports that this chamber has Euler characteristic 1. The other two chambers each have one hole and thus have Euler characteristic 0. ","category":"page"},{"location":"#Example:-3-random-polynomials","page":"Chambers.jl","title":"Example: 3 random polynomials","text":"","category":"section"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"We can set up a random example as follows.","category":"page"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"julia> using Chambers\njulia> @var v[1:3];\njulia> f = [rand_poly(Float64, v, d) for d in [2, 3, 3]];\njulia> C = chambers(f)","category":"page"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"Here, f consists of 3 random polynomials in v. The degrees of these polynomials are [2, 3, 3]. The coefficients are chosen from a Gaussian distribution.","category":"page"},{"location":"#Documentation:-Output","page":"Chambers.jl","title":"Documentation: Output","text":"","category":"section"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"ChambersResult\nChamber","category":"page"},{"location":"#Chambers.ChambersResult","page":"Chambers.jl","title":"Chambers.ChambersResult","text":"ChambersResult\n\nA struct that collects all chambers in the complement of a hypersurface arragement.\n\n\n\n\n\n","category":"type"},{"location":"#Chambers.Chamber","page":"Chambers.jl","title":"Chambers.Chamber","text":"Chamber\n\nA struct that contains all information about a chamber.\n\n\n\n\n\n","category":"type"},{"location":"#Documentation:-Main-functions","page":"Chambers.jl","title":"Documentation: Main functions","text":"","category":"section"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"chambers\naffine_chambers\nprojective_chambers","category":"page"},{"location":"#Chambers.chambers","page":"Chambers.jl","title":"Chambers.chambers","text":"chambers(C::ChambersResult)\n\nReturns the vector of chambers in C.\n\n\n\n\n\nchambers(f::Vector{Expression})\nchambers(f::System)\n\nInput a list of hypersurfaces 'f = [f1,...fk]'. Outputs the chambers in the complement of the hypersurface arrangement, whether they are bounded or not, their sign patterns, Euler characteristic and the indices of the critical points in each chamber.\n\nOptions:\n\nshow_progress = true: if true, prints the progress of the computation to the terminal.\nprojective_fusion = true: if true, the algorithm computes which of the chambers are fused at infinity.\ns: a list of integers [s_1, ..., s_k, s_{k+1}] such that s_1, ..., s_k>0, s_{k+1}<0 and 2 s_{k+1} > s_1 deg(f_1) + ... + s_k deg(f_k).\nepsilon = 1e-6: how close from each critical point do we do the path tracking.\nreltol = 1e-6, abstol = 1e-9: parameters for the accuracy of the ODE solver.\nmonodromy_options = MonodromyOptions(max_loops_no_progress = 25): pass options for monodromy.\nstart_pair_using_newton::Bool = false: if true, the algorithm tries to compute a start pair for monodromy by using Newton's methods. Can reduce the number of critical points, but is less stable.\n\nExample\n\nusing Chambers\n@var x y\nf = [x^2 + y^2 - 1; x^2 + y^2 - 4];\nchambers(f)\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.affine_chambers","page":"Chambers.jl","title":"Chambers.affine_chambers","text":"affine_chambers(f::Vector{Expression})\naffine_chambers(f::System)\n\nInput a list of affine hypersurfaces 'f = [f1,...fk]'.  Outputs the chambers in the complement of the hypersurface arrangement, their sign patterns, Euler characteristic and the indices of the critical points in each chamber. Accepts the same options as chambers.\n\nExample\n\nusing Chambers\n@var x y\nf = [x^2 + y^2 - 1; x^2 + y^2 - 4];\naffine_chambers(f)\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.projective_chambers","page":"Chambers.jl","title":"Chambers.projective_chambers","text":"projective_chambers(C::ChambersResult)\n\nReturns a vector of vectors. The entries of the vectors correspond to those chambers in C, which are fused in projective space.\n\nThe following code will return a vector of vectors of type Chamber:\n\nusing Chambers\n@var x y\nf = [x^2 + y^2 - 1; x^2 + y^2 - 4];\nC = chambers(f)\np = projective_chambers(C)\nmap(pᵢ -> chambers(C)[pᵢ], p)\n\n\n\n\n\n","category":"function"},{"location":"#Documentation:-Helper-functions","page":"Chambers.jl","title":"Documentation: Helper functions","text":"","category":"section"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"Functions to call on Chamber.","category":"page"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"χ\nμ\ncritical_points\nis_bounded","category":"page"},{"location":"#Chambers.χ","page":"Chambers.jl","title":"Chambers.χ","text":"χ(C::Chamber)\n\nReturns the Euler characteristic.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.μ","page":"Chambers.jl","title":"Chambers.μ","text":"μ(C::Chamber)\n\nReturns the index vector.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.critical_points","page":"Chambers.jl","title":"Chambers.critical_points","text":"critical_points(C::Chamber)\n\nReturns the critical points in C.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.is_bounded","page":"Chambers.jl","title":"Chambers.is_bounded","text":"is_bounded(C::Chamber)\n\nReturns a boolean that is true, if C is bounded. If that information was not computed, it simply returns nothing.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"Functions to call on ChambersResult.","category":"page"},{"location":"","page":"Chambers.jl","title":"Chambers.jl","text":"nchambers\nnbounded\nnunbounded\nbounded\nunbounded\neuler_characteristics\nindex_vectors\nncritical_complex\nncritical_real\ng","category":"page"},{"location":"#Chambers.nchambers","page":"Chambers.jl","title":"Chambers.nchambers","text":"nchambers(C::ChambersResult)\n\nReturns the number of chambers in F.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.nbounded","page":"Chambers.jl","title":"Chambers.nbounded","text":"nbounded(C::ChambersResult)\n\nReturns the number of bounded chambers in C.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.nunbounded","page":"Chambers.jl","title":"Chambers.nunbounded","text":"nunbounded(C::ChambersResult)\n\nReturns the number of unbounded chambers in C.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.bounded","page":"Chambers.jl","title":"Chambers.bounded","text":"bounded(C::ChambersResult)\n\nReturns the bounded chambers in C.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.unbounded","page":"Chambers.jl","title":"Chambers.unbounded","text":"unbounded(C::ChambersResult)\n\nReturns the unbounded chambers in C.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.euler_characteristics","page":"Chambers.jl","title":"Chambers.euler_characteristics","text":"euler_characteristics(C::ChambersResult)\n\nReturns the Euler characteristics of the chambers in C.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.index_vectors","page":"Chambers.jl","title":"Chambers.index_vectors","text":"euler_characteristics(C::ChambersResult)\n\nReturns the index vectors of the chambers in C.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.ncritical_complex","page":"Chambers.jl","title":"Chambers.ncritical_complex","text":"ncritical_complex(C::ChambersResult)\n\nReturns the number of complex critical points of the Morse function in C.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.ncritical_real","page":"Chambers.jl","title":"Chambers.ncritical_real","text":"ncritical_real(C::ChambersResult)\n\nReturns the number of real critical points of the Morse function in C.\n\n\n\n\n\n","category":"function"},{"location":"#Chambers.g","page":"Chambers.jl","title":"Chambers.g","text":"g(C::ChambersResult)\n\nReturns the Morse function in C.\n\n\n\n\n\n","category":"function"}]
}