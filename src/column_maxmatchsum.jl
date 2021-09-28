struct one_matching_output
    q::Float64
    mi::Vector{Int64}
end

struct all_matching_output
    q::Vector{Float64}
    mi::Vector{Int64}
    mj::Vector{Int64}
    medges::Int64
end

function intmatch(n::Int64,m::Int64,nedges::Int64,
                v1::Vector{Int64},v2::Vector{Int64},weight::Vector{Float64})


l1 = zeros(Float64,n)
l2 = zeros(Float64,n+m)
s = ones(Int64,n+m)
t = -1*ones(Int64,n+m)
offset = zeros(Int64,n)
deg = ones(Int64,n)
list = zeros(Int64,nedges+n)
index = zeros(Int64,nedges+n)
w = zeros(Float64,nedges+n)
match1 = -1*ones(Int64,n)
match2 = -1*ones(Int64,n+m)
tmod = zeros(Int64,m+n)
ntmod = 0

for i = 1:nedges
    deg[v1[i]] += 1
end

for i  = 2:n
    offset[i] = offset[i-1] + deg[i-1]
end
offset .+= 1

deg .= 0

for i  = 1:nedges
    list[offset[v1[i]] + deg[v1[i]]] = v2[i]
    w[offset[v1[i]] + deg[v1[i]]] = weight[i]
    index[offset[v1[i]] + deg[v1[i]]] = i
    deg[v1[i]] += 1
end

for i = 1:n
    list[offset[i]+deg[i]] = m + i
    w[offset[i]+deg[i]] = 0
    index[offset[i]+deg[i]] = -1
    deg[i] += 1
end

for i = 1:n
    for j = 0:deg[i]-1
        if w[offset[i]+j] > l1[i]
            l1[i] = w[offset[i]+j]
        end
    end
end

i = 1
while i<=n
    for j = 1:ntmod
        t[tmod[j]] = -1
    end
    ntmod = 0
    p = 1
    q = 1
    s[1] = i
    while p<=q
        if match1[i] >= 1
            break
        end
        k = s[p]
        for r = 0:deg[k]-1
            j = list[offset[k]+r]
            val1totest = w[offset[k]+r]
            val2totest = l1[k] + l2[j] - 1e-8
            if w[offset[k]+r] < val2totest
                continue
            end
            if t[j] < 0
                q = q + 1
                s[q] = match2[j]
                t[j] = k
                ntmod += 1
                tmod[ntmod] = j
                if match2[j] < 0
                    while j>=1
                        k = match2[j] = t[j]
                        p = match1[k]
                        match1[k] = j
                        j = p
                    end
                    break
                end
            end
        end
        p += 1
    end
    p -= 1

    if match1[i] < 0
        al = 1e20
        for j = 1:p
            t1 = s[j]
            for k = 0:deg[t1]-1
                t2 = list[offset[t1] + k]
                if t[t2] < 0 && l1[t1] + l2[t2] - w[offset[t1] + k] < al
                    al = l1[t1] + l2[t2] - w[offset[t1] + k]
                end
            end
        end
        for j = 1:p
            vtemp = s[j]
            vvtemp = l1[s[j]]
            l1[s[j]] -= al
            vvtemp = l1[s[j]]
        end
        for j = 1:ntmod
            l2[tmod[j]] += al
        end
    else
        i += 1
    end
end


ret = 0.0
for i = 1:n
    for j = 0:deg[i] - 1
        if list[offset[i] + j] == match1[i]
            ret += w[offset[i]+j]
        end
    end
end

mi = zeros(Int64,nedges)
for i = 1:n
    if match1[i]<=m
        for j = 0:deg[i]-1
            if list[offset[i] + j] == match1[i]
                mi[index[offset[i]+j]] = 1
            end
        end
    end
end

return one_matching_output(ret,mi)
end


function column_maxmatchsum(M::Int64,N::Int64,Qp::Vector{Int64},Qr::Vector{Int64},Qv::Vector{Float64},
                                m::Int64,n::Int64,nedges::Int64,li::Vector{Int64},lj::Vector{Int64})

minnm = n
if m<minnm
    minnm = m
end

q = zeros(Float64,N)
mi = zeros(Int64,Qp[N]-1)
mj = zeros(Int64,Qp[N]-1)

medges = 1

lwork1 = ones(Int64,m)
lind1 = -1*ones(Int64,m)
lwork2 = ones(Int64,n)
lind2 = -1*ones(Int64,n)

max_col_nonzeros = 0
for j = 1:N
    col_nonzeros = Qp[j+1] - Qp[j]
    if col_nonzeros >= max_col_nonzeros
        max_col_nonzeros = col_nonzeros
    end
end

se1 = zeros(Int64,max_col_nonzeros)
se2 = zeros(Int64,max_col_nonzeros)
sw = zeros(Float64,max_col_nonzeros)
sqi = zeros(Int64,max_col_nonzeros)
sqj = zeros(Int64,max_col_nonzeros)

for j = 1:N
    smalledges = 1
    nsmall1 = 1
    nsmall2 = 1

    for nzi = Qp[j]:Qp[j+1]-1
        i = Qr[nzi]
        v1 = li[i]
        v2 = lj[i]
        sv1 = -1
        sv2 = -1
        if lind1[v1] < 0
            sv1 = nsmall1
            lind1[v1] = sv1
            lwork1[sv1] = v1
            nsmall1 += 1
        else
            sv1 = lind1[v1]
        end
        if lind2[v2] < 0
            sv2 = nsmall2
            lind2[v2] = sv2
            lwork2[sv2] = v2
            nsmall2 += 1
        else
            sv2 = lind2[v2]
        end
        se1[smalledges] = sv1
        se2[smalledges] = sv2
        sw[smalledges] = Qv[nzi]
        sqi[smalledges] = i
        sqj[smalledges] = j
        smalledges += 1
    end
    nsmall1 -= 1
    nsmall2 -= 1
    smalledges -= 1

    if smalledges==0
        q[j] = 0
        continue
    end

    one_matching = intmatch(nsmall1,nsmall2,smalledges,se1,se2,sw)
    q[j] = one_matching.q
    smi = one_matching.mi

    for k = 1:smalledges
        if smi[k] > 0
            mi[medges] = sqi[k]
            mj[medges] = sqj[k]
            medges += 1
        end
    end

    for k = 1:nsmall1
        lind1[lwork1[k]] = -1
    end

    for k = 1:nsmall2
        lind2[lwork2[k]] = -1
    end

end

medges -= 1
return all_matching_output(q,mi,mj,medges)

end
