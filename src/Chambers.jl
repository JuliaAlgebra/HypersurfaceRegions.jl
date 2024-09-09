module Regions

import LinearAlgebra, LightGraphs, Random
import HomotopyContinuation
import DifferentialEquations, SciMLBase
import ProgressMeter
using PrettyTables, Crayons

const HC = HomotopyContinuation
const DE = DifferentialEquations
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
include("regions_main.jl")

end
