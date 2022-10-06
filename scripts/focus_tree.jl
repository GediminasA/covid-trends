#!/usr/bin/env julia
using ArgMacros
using NewickTree


function get_labels_from_file(nodesfile::String)::Array{String,1}
    labels_text = read(nodesfile, String)
    labelstmp = split(labels_text)
    labels = Array{String,1}()
    labelstmp = filter((i) -> i != " ", labelstmp)
    labelstmp = filter((i) -> i != "", labelstmp)
    removechar = ['\n', '\r', ' ', '\t']
    for l in labelstmp
        cleanl = replace(l,removechar => "")
        push!(labels,cleanl)
    end
    return(labels)
end

function make_focused!(t::NewickTree.Node,focus::Array{String,1};leave_siblings = true)
    focused_lvs_names = focus
    not_focused_lvs = Array{Node{UInt16, NewickData{Float64, String}},1}()
    focused_lvs = Array{Node{UInt16, NewickData{Float64, String}},1}()

    for l in NewickTree.getleaves(t)
        if !(l.data.name in focused_lvs_names)
            push!(not_focused_lvs,l)
        else
            push!(focused_lvs,l)
        end
    end
    nodes_4_deletion = Array{Node{UInt16, NewickData{Float64, String}},1}()
    for n in postwalk(t)
        l_lvs = getleaves(n)
        #println(n.data.name)
        #println(l_lvs)
        #println("___")
        cmn = length(intersect(l_lvs,focused_lvs))
        if cmn == 0 && (!NewickTree.isleaf(n) || !leave_siblings) && !NewickTree.isroot(n)
            push!(nodes_4_deletion,n)
        end
    end
    for n in nodes_4_deletion
        p = n.parent
        NewickTree.delete!(p,n)
    end

    for l in t.children
        if l in not_focused_lvs
            NewickTree.delete!(t,l)
        end
    end
end

function main()
    @inlinearguments begin
        @argumentrequired String intree  "-i" "--innwk"
        @arghelp "Input tree in newick tree"
        @argumentrequired String outtree "-o" "--onnwk"
        @arghelp "Output tree in newick tree"
        @argumentflag keepsiblings "-s" "--keep-sibling"
        @arghelp "Keep sibling"
        @argumentrequired String nodesfile "-l" "Labels for focus" 
        @arghelp "List of leaf names to focus on"
    end
    tini = readnw(readline(intree))
    labels4focus = get_labels_from_file(nodesfile)
    make_focused!(tini,labels4focus,leave_siblings = keepsiblings)
    FF = open(outtree,"w")
    tr = nwstr(tini,internal=true)
    print(FF,tr)
    close(FF)
end
main()
