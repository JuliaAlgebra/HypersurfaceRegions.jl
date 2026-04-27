module HypersurfaceRegions

import LinearAlgebra, LightGraphs, Random
import HomotopyContinuation
import OrdinaryDiffEq, SciMLBase
import ProgressMeter
using PrettyTables, Crayons

const HC = HomotopyContinuation
const DE = OrdinaryDiffEq
const LA = LinearAlgebra
const LG = LightGraphs
const PM = ProgressMeter

using Reexport: @reexport
@reexport using HomotopyContinuation


include("progressmeter.jl")
include("output.jl")
include("hessian.jl")
include("path_tracking.jl")
include("partition.jl")
include("affine_regions.jl")
include("membership.jl")
include("regions.jl")

end
