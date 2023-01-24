#!/usr/bin/env julia
using ArgMacros
using NewickTree
using DataFrames
using Query
using SHA
using JSON

function parse_country_data(injson::String)::Dict{String,String}
    node_country = Dict{String,String}()
    fj = open(injson,"r")
    strj = read(fj,String)
    js = JSON.parse(strj)["nodes"]
    for n in keys(js)
        node_country[n] = String(js[n]["country"])
    end
    return node_country
end
function get_country(name::String,node_country::Dict{String,String})
    if name in keys(node_country)
        return(node_country[name])
    else 
        return("NA")
    end
end


function get_cons_pair2(t::NewickTree.Node,get_country::Function,node_country::Dict{String,String})::String
    out = "NO"
    for n in postwalk(t)
        internal = true
        if !NewickTree.isleaf(n) 
            for l in n.children
                if NewickTree.isleaf(l)
                    internal = false
                    break
                end
            end
        else
            internal = false
        end
        if internal &&length(n.children) == 1 && !NewickTree.isroot(n) && ! (n.parent.data.name == "")
             ch = n.children[1]
             n_ch_cntr = get_country(ch.data.name,node_country)
             n_p_cntr = get_country(n.data.name,node_country)
             if n_ch_cntr == n_p_cntr
                 out = n.data.name
             end
        end
    end
    return(out)
end

function remove_2degree_node(t::NewickTree.Node,nodename::String,node_country::Dict{String,String})
    out = deepcopy(t)
    tonode = NaN
    found = false
    for p in postwalk(out)
        if p.data.name == nodename
            found = true
            if length(p.children) >1 
                println("ERR...")
                throw(UndefVarError(:nodename))
            else
                println("Removinggg, ",nodename)
#                 println(p.data.name,"  ",nodename," ", p.data.name == nodename)
                ch = deepcopy(p.children[1])
                ch.data.distance = ch.data.distance + p.data.distance
                p0 = p.parent
                println("with parent, ",p0.data)
                tonode = ch.data.name
                ch.parent = p0
#                 println("UUUUUUUUUU0")
#                 for c in p0.children 
#                     println(c.data.name)
#                 end
#                 println("UUUUUUUUUU1")
                NewickTree.delete!(p0,p)
                push!(p0.children,ch)
#                 println("UUUUUUUUUU1")
#                 for c in p0.children 
#                     println(c.data.name)
#                 end
#                 println("UUUUUUUUUU2")
               
            end
        end
    end
#     println("recheck")
#     for p in postwalk(out)
#         println(p.data.name)
#         if p.data.name == nodename || p.data.name in ["NODE_0000051","NODE_0000034"]
#             println("HHHHHHHHHHHHHHHHHHHH")
#             println("**",p.parent.data.name)
#             for c in p.children
#                 println("- ",c.data.name)
#             end
#         end
        
#     end
#     #println("&&&&&&&")
    if !found
        println(nodename, " not found...exiting")
        throw(UndefVarError(:nodename))
    end
    return(out,tonode)
end

function remove_2dg(tt::NewickTree.Node,get_country::Function,node_country::Dict{String,String}) 
        data_joined = Dict{String,Array{String,1}}()
        ttl = deepcopy(tt)
        n4del = get_cons_pair2(ttl,get_country,node_country)
        while n4del != "NO"
            ttl,tonode = remove_2degree_node(ttl,n4del,node_country)
            #println([n4del,tonode])
            for nn in [n4del,tonode]
                if !(nn in keys(data_joined))
                    data_joined[nn] = Array{String,1}([nn])
                end
            end
            data_joined[tonode] = vcat(data_joined[tonode],data_joined[n4del])
            
            #collecdata on merge
            n4del = get_cons_pair2(ttl,get_country,node_country)
            #collecdata on merge
        end
        return(ttl,data_joined)
end




function main()
    @inlinearguments begin
        @argumentrequired String intree  "-i" "--innwk"
        @arghelp "Input tree in newick tree"
        @argumentrequired String indata  "-j" "--indata"
        @arghelp "Input tree data"
        @argumentrequired String outtree "-o" "--outnwk"
        @arghelp "Output tree"
        @argumentrequired String outdata "-d" "--outdata"
        @arghelp "Output country data"
    end
    # read in tree 
    node_country=parse_country_data(indata)
    tree = readnw(readline(intree))
    ttout , data_joined= remove_2dg(tree,get_country,node_country)    
    dkeys = keys(data_joined)
    for n in postwalk(ttout)
        if NewickTree.isroot(n)
            n.data.name = "ROOT"
        end
        nn = n.data.name
        if !(nn in dkeys)
            #println("A>|",nn)
            data_joined[nn] = [nn]
        end
    end
    
    FF=open(outtree,"w")
    tr = nwstr(ttout,internal=true)
    print(FF,tr)
    close(FF)
    F=outdata 
    FF = open(F,"w")
    print(data_joined)
    JSON.print(FF, data_joined)
    close(FF)




    # for k in keys(data_joined)
    #     println(k," ",length(data_joined[k])," ",data_joined[k])
    # end
    #println(get_country("NODE_0000833",node_country))
end
main()
# mamba install -c anaconda cython python=3.8 ; mamba install gurobi