module NetworkAlign

using MatrixNetworks
using SparseArrays

import Printf: @printf
import DelimitedFiles: readdlm
import LinearAlgebra: triu, dot, norm, issymmetric
import KahanSummation: sum_kbn

"""
Module ``NetworkAlign``: Documentation on the module

You can check the readme file here: \n
"https://github.com/nassarhuda/NetworkAlign.jl/blob/master/README.md"
"""
#NetworkAlign

include("normout.jl")
include("othermaxplus.jl")
include("othersum.jl")
include("round_messages.jl")
include("column_maxmatchsum.jl")
include("load_netalign_problem.jl")
include("netalignmr.jl")
include("isorank.jl")
include("netalignbp.jl")
include("count_overlap.jl")
include("make_squares.jl")
include("netalign_setup.jl")

export load_netalign_problem, netalignmr, isorank, netalignbp,
        count_overlap, make_squares, netalign_setup, build_Si

end # end module
