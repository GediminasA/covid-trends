#!/usr/bin/env julia
using ArgMacros
using NewickTree
using DataFrames
using Query


function main()
    @inlinearguments begin
        @argumentrequired String intree  "-i" "--innwk"
        @arghelp "Input tree in newick tree"
        @argumentrequired String outgroup "-o" "--od"
        @arghelp "file to extract outgroup for reroooting other tree"
    end
    # read in tree 
    tini = readnw(readline(intree))
    parts = Array{Array{String,1},1}()
    for n in postwalk(tini)
        if NewickTree.isroot(n)
            for l in n.children
                l_kids = getleaves(l)
                lnsms = Array{String,1}()
                for lk in l_kids
                    push!(lnsms,lk.data.name)
                end
                push!(parts,lnsms)
            end
        end
    end
    sort!(parts,by= x ->  length(x))
    fo=open(outgroup,"w")
    for n in parts[1]
        println(fo,n)
    end
    close(fo)
end
main()
# mamba install -c anaconda cython python=3.8 ; mamba install gurobi