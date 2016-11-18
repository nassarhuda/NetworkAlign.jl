module NetworkAlign

using MatrixNetworks

"""
Module ``NetworkAlign``: Documentation on the module

You can check the readme file here: \n
"https://github.com/nassarhuda/NetworkAlign.jl/blob/master/README.md"
"""
NetworkAlign

include("normout.jl")
include("othermaxplus.jl")
include("othersum.jl")
include("round_messages.jl")
include("column_maxmatchsum.jl")
include("load_netalign_problem.jl")
include("netalignmr.jl")
include("isorank.jl")
include("netalignbp.jl")

export load_netalign_problem, netalignmr, isorank, netalignbp

end # end module
