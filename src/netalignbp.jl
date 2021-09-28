"""
NETALIGNBP
----------
    solve the network alignment problem with Belief Propagation
"""

function netalignbp(S::SparseMatrixCSC{T,Int64},w::Vector{F},
                            a::R,b::K,li::Vector{Int64},
                            lj::Vector{Int64},gam::Float64,
                            dtype::Int64,maxiter::Int64,verbose::Bool) where {T,F,R,K}
  nedges = length(li)
  nsquares = nnz(S)/2
  m = maximum(li)
  n = maximum(lj)
  SU = findall(x->x!=0,triu(S))#only consider each square once.
  
  #compute a vector that allows us to transpose data between squares.
  #Recall the BP algorithm requires edge(i,j) to send messages to edge(r,s)
  #when (i,j) and (r,s) form a square.  However, edge(i,j) needs information
  #the information that (r,s) sent it from the previous iteraton.
  #If we imagine that each non-zero of S has a value, we want to tranpose
  #these values.

  SI = spzeros(Int,size(S,1),size(S,1))
  SI[SU] = 1:length(SU)
  # this is faster than the old version of ind2sub, now SU is an array of CartesianIndex type
  # SI = sparse(sui,suj,1:length(sui),size(S,1),size(S,2)) #assign indices
  
  SI = SI + SI' #SI now has symmetric indices
  si = SI.rowval
  sind = SI.nzval

  # faster than:
  # (si, sj) = ind2sub(size(SI),find(SI))
  # sind = nonzeros(SI)
  
  SP = sparse(si,sind,true,size(S,1),nsquares); # each column in SP has 2 nz
  
  sij_cartesian = findall(x->x!=0,SP)
  sij = [x.I[1] for x in sij_cartesian]
  sijrs = [x.I[2] for x in sij_cartesian]


  spair = sind
  len_spair = length(spair)
  spair[1:2:end] = 2:2:len_spair
  spair[2:2:end] = 1:2:len_spair

  # sij, sijrs maps between rows and squares
  # spair is now an indexing vector that accomplishes what we need

  # Initialize the messages
  ma = zeros(Float64,nedges)
  mb = zeros(Float64,nedges)
  ms = zeros(Float64,nnz(S))
  sums = zeros(Float64,nedges)
  damping = gam
  curdamp = 1
  iter = 1
  alpha = a
  beta = b

  # Initialize history
  hista = zeros(Float64,maxiter,4) # history of messages from ei->a vertices
  histb = zeros(Float64,maxiter,4) # history of messages from ei->b vertices
  fbest = 0
  fbestiter = 0
  if verbose # print the header
    @printf("%4s   %4s   %7s %7s %7s %7s   %7s %7s %7s %7s\n",
        "best", "iter", "obj_ma", "wght_ma", "card_ma", "over_ma",
            "obj_mb", "wght_mb", "card_mb", "over_mb")
  end
  # setup the matching problem once
  Matching_Setup = bipartite_matching_setup(w,li,lj,m,n)
  tripi = Matching_Setup.tripi
  matm = Matching_Setup.m
  matn = Matching_Setup.n
  rp = Matching_Setup.rp
  ci = Matching_Setup.ci
  mperm = tripi[tripi.>0]

  while iter <= maxiter
    prevma = copy(ma)
    prevmb = copy(mb)
    prevms = copy(ms)
    prevsums = copy(sums)
    curdamp = damping * curdamp

    omaxb = max.((othermaxplus(2,li,lj,mb,m,n)),0)
    omaxa = max.((othermaxplus(1,li,lj,ma,m,n)),0)

    msflip = ms[spair] # swap ij->ijrs to rs->ijrs

    mymsflip = msflip .+ beta
    mymsflip = min.(beta,mymsflip)
    mymsflip = max.(0,mymsflip)

    sums = zeros(Float64,nedges)
    for i = 1:length(sij)
        sums[sij[i]] += mymsflip[i]
    end

    ma = alpha*w .- omaxb .+ sums
    mb = alpha*w .- omaxa .+ sums

    ms = alpha*w[sij]-(omaxb[sij] + omaxa[sij])
    ms = ms .+ othersum(sij,mymsflip,nedges)
    # Original updates
    if dtype == 1
      ma = curdamp*(ma) .+ (1-curdamp)*(prevma)
      mb = curdamp*(mb) .+ (1-curdamp)*(prevmb)
      ms = curdamp*(ms) .+ (1-curdamp)*(prevms)
    elseif dtype == 2
      ma = ma .+ (1-curdamp)*(prevma.+prevmb.-alpha*w+prevsums)
      mb = mb .+ (1-curdamp)*(prevmb.+prevma.-alpha*w+prevsums)
      ms = ms .+ (1-curdamp)*(prevms.+prevms[spair].-beta)
    elseif dtype == 3
      ma = curdamp*ma .+ (1-curdamp)*(prevma.+prevmb.-alpha*w.+prevsums)
      mb = curdamp*mb .+ (1-curdamp)*(prevmb.+prevma.-alpha*w.+prevsums)
      ms = curdamp*ms .+ (1-curdamp)*(prevms.+prevms[spair].-beta)
    end
    
    # now compute the matchings
    hista[iter,:] = round_messages(ma,S,w,alpha,beta,rp,ci,tripi,matn,matm,mperm,li,lj)
    histb[iter,:] = round_messages(mb,S,w,alpha,beta,rp,ci,tripi,matn,matm,mperm,li,lj)
    if hista[iter,1] > fbest
      fbestiter = iter
      mbest = ma
      fbest = hista[iter,1]
    end
    if histb[iter,1] > fbest
      fbestiter = -iter
      mbest = mb
      fbest = histb[iter,1]
    end
    if verbose
      if fbestiter == iter
        bestchar = "*a"
      elseif fbestiter == -iter
        bestchar = "*b"
      else
        bestchar = ""
      end
      @printf("%4s   %4i    %5.2f %5.2f %4d %7d %14.2f %5.2f %4d %7d\n",
            bestchar, iter, hista[iter,1], hista[iter,2], hista[iter,3], hista[iter,4],
              histb[iter,1], histb[iter,2], histb[iter,3], histb[iter,4]);

    end
    iter += 1
  end
  return (ma,S,w,alpha,beta,rp,ci,tripi,matn,matm,mperm,li,lj)

end
