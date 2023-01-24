#!/usr/bin/env julia
using ArgMacros
using NewickTree
using DataFrames
using Query
using SHA

function main()
    @inlinearguments begin
        @argumentrequired String intree  "-i" "--innwk"
        @arghelp "Input tree in newick tree"
        @argumentrequired String outdir "-o" "--od"
        @arghelp "Output directory"
        @argumentdefault Int64 4000 maximumnumber  "-m" "--max"
        @arghelp "Maximum number"
        @argumentdefault Int64 10 minimumnumber "-n" "--min"
        @arghelp "Minimum number"
    end
    # read in tree 
    tini_name = intree
    maxn = maximumnumber
    minn = minimumnumber
    exhausted = false
    mkpath(outdir)
    ct = 0
    tini = readnw(readline(tini_name))
    covered_nodes = Array{String,1}()

    cnt = 0
    # add names to nodes just in case there are no
    for n in postwalk(tini)
        if !NewickTree.isleaf(n) && !NewickTree.isroot(n) 
            l_lvs = getleaves(n)
            names = Array{String,1}()
            for l in l_lvs
                push!(names,l.data.name)
            end
            n.data.name = bytes2hex(sha256(join(names)))
        end
    end
    while !exhausted
        cnt += 1
        data = DataFrame(node=Array{String,1}(),size=Array{Int64,1}())
        for n in postwalk(tini)
            if !NewickTree.isleaf(n) && !NewickTree.isroot(n) 
                l_lvs = getleaves(n)
                push!(data,[n.data.name,length(l_lvs)])
            end
        end
        data = @from i in data begin
            @orderby descending(i.size)
            @select i
            @collect DataFrame
        end
        #println(first(data,10))
        data2 = @from i in data begin
            @orderby descending(i.size)
            @where i.size <= maxn 
            @select i
            @collect DataFrame
        end
        #println(first(data2,10))
        data3 = @from i in data2 begin
            @orderby descending(i.size)
            @where !(i.node in covered_nodes)
            @select i
            @collect DataFrame
        end
        
        if nrow(data3) < 3
            println("No tips left after filtring ...exiting")
            println(data3)
            break
        end
        chosenid = data3.node[1]
        node = NA
        ct = ct + 1
        for n in postwalk(tini)
            if !NewickTree.isleaf(n) && !NewickTree.isroot(n) && n.data.name == chosenid
                node=n
            end
        end
        
        #get all leaves
        for n in postwalk(node)
            push!(covered_nodes,n.data.name)
        end
        
        nl = length(getleaves(node))
        
        if nl >  minn
            outname = outdir*"/" * "tree_" * "$ct"*".nwk"
            F = outname
            FF = open(F,"w")
            tr = nwstr(node,internal=false)
            print(FF,tr)
            close(FF)
            println(stderr," leaves number is $nl , writingf out to $outname")
        else
            exhausted = true
            print(" leaves number is $nl stopping...")
        end
    end
end
main()
# mamba install -c anaconda cython python=3.8 ; mamba install gurobi