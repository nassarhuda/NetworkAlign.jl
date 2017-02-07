# Netalign
Network Alignment Algorithms in Julia
* sample runs:
```
Pkg.clone("https://github.com/nassarhuda/NetworkAlign.jl")
using NetworkAlign
```

### dataset name: example-overlap 
### size of the problem: 12 nodes
```
S,w,li,lj,A,B,L = load_netalign_problem("example-overlap")
a = 1;
b = 1;
stepm = 25;
rtype = 1;
maxiter = 10;
verbose = true;
gamma = 0.4;
xbest,st,status,hist = netalignmr(S,w,a,b,li,lj,gamma,stepm,rtype,maxiter,verbose)
x,flag,reshist = isorank(S,w,a,b,li,lj)
x,flag,reshist = netalignbp(S,w,a,b,li,lj,0.99,2,100,true)
```

### dataset name: lcsh2wiki-small
### size of the problem: 16952 nodes
```
S,w,li,lj,A,B,L = load_netalign_problem("lcsh2wiki-small")
a = 1;
b = 1;
stepm = 25;
rtype = 1;
maxiter = 100;
verbose = true;
gamma = 0.4;
xbest,st,status,hist = netalignmr(S,w,a,b,li,lj,gamma,stepm,rtype,maxiter,verbose)
x,flag,reshist = isorank(S,w,a,b,li,lj)
x,flag,reshist = netalignbp(S,w,a,b,li,lj,0.99,2,100,true)
```

### dataset name: lcsh2wiki-full
### size of the problem: 4971629 nodes
```
S,w,li,lj,A,B,L = load_netalign_problem("lcsh2wiki-full")
a = 1;
b = 1;
stepm = 25;
rtype = 1;
maxiter = 100;
verbose = true;
gamma = 0.4;
xbest,st,status,hist = netalignmr(S,w,a,b,li,lj,gamma,stepm,rtype,maxiter,verbose)
x,flag,reshist = isorank(S,w,a,b,li,lj)
x,flag,reshist = netalignbp(S,w,a,b,li,lj,0.99,2,100,true)
```
## we can start with networks A and B, and similarity matrix L, and construct the alignment problem as follows:
### dataset name: example-2
### size of the problem: 13 nodes
```
# load data
mainloc = joinpath(Pkg.dir("NetworkAlign"),"data/example-2")
location = join([mainloc,"_A.smat"])
A = MatrixNetworks.readSMAT(location)
location = join([mainloc,"_B.smat"])
B = MatrixNetworks.readSMAT(location)
location = join([mainloc,"_L.smat"])
L = MatrixNetworks.readSMAT(location)

# setup experiment
S,w,li,lj = netalign_setup(A,B,L)

# align networks
a = 1;
b = 1;
stepm = 25;
rtype = 1;
maxiter = 10;
verbose = true;
gamma = 0.4;
xbest,st,status,hist = netalignmr(S,w,a,b,li,lj,gamma,stepm,rtype,maxiter,verbose)
x,flag,reshist = isorank(S,w,a,b,li,lj)
x,flag,reshist = netalignbp(S,w,a,b,li,lj,0.99,2,100,true)
```
