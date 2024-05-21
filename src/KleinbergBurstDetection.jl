module KleinbergBurstDetection

using AbstractTrees
using Printf

include("model_costs.jl")
include("dynamic_program.jl")

export BurstNode, BurstStructure, detect_bursts, hierarchical_burst_structure

mutable struct BurstNode
    activity::Int64
    start_index::Int64
    end_index::Int64
    start_time::Float64
    end_time::Float64
    children::Vector{BurstNode}
    parent::Union{Nothing,BurstNode}
end

AbstractTrees.children(x::BurstNode) = x.children
AbstractTrees.nodevalue(x::BurstNode) = x.activity

AbstractTrees.printnode(io::IO,x::BurstNode) = @printf(io,"Lvl: %i start: %.5f end: %.5f",x.activity,x.start_time,x.end_time)

# printing the tree structure
function Base.show(io::IO,x::BurstNode)
    print_tree(io,x)
end

struct BurstStructure
    times::Vector{Float64}
    activity::Vector{Int64}
    hierarchy::BurstNode
end

function Base.show(io::IO,x::BurstStructure)
    print(io,"KleinbergBurstDetection.BurstStructure\n")

    print(io,"events: ")
    print(io,x.times)
    print(io,"\n")

    print(io,"activity state: ")
    print(io,x.activity)
    print(io,"\n")

    print(io,"hierarchy: \n")
    print(io,x.hierarchy)
end

####
# Given the output of Kleinbergs algorithm this constructs an AbstractTrees.jl tree representing the hierarchical structure of bursts.
####
function hierarchical_burst_structure(times::Vector{Float64},state_sequence::Vector{Int64})
    n = size(state_sequence,1)
    # create the root node
    root = BurstNode(0,1,n,times[begin],times[end],[],nothing)
    # now iterate through the states and build the tree
    cursor = root
    for t ∈ 1:n
        while cursor.activity != state_sequence[t]
            if cursor.activity < state_sequence[t]
                # make a new intermediate node
                next_cursor = BurstNode(cursor.activity+1,t,t,times[t],times[t],[],cursor)
                push!(cursor.children,next_cursor)
                cursor = next_cursor
            else
                next_cursor = cursor.parent
                cursor.end_index = t
                cursor.end_time = times[t]
                cursor = next_cursor
            end
        end
        if cursor.end_index < t
            cursor.end_index = t
            cursor.end_time = times[t]
        end
    end
    return root
end

###
# printing the tree structure
###

###
# Given the list of event times this function performs Kleinbergs Algorithm and 
# associates each event time with an activity level where 0 is no burst and >=1 is some burst
# activity is increasing exponentially with the activity level with base s.
# The parameter γ describes how many bursts are detected, γ=1.0 is a rather unsensitive filter. If one
# uses 0 < γ < 1 the more burst events are detected the smaller γ is. 
###
function detect_bursts(event_times::Vector{Float64},s=2.0,γ=1.0)

    # check input values
    if s <= 1.0
        throw(ArgumentError("[KBD ERROR] the scaling factor s needs to be larger than 1.0."))
    end
    if size(event_times,1) < 3
        throw(ArgumentError("[KBD ERROR] there need to be at least 3 events in the time series."))
    end

    # make sure events are sorted
    order = sortperm(event_times)
    sorted_events = event_times[order]

    # rescale event times to [0,1]
    n = size(sorted_events,1)
    min_t = minimum(sorted_events)
    scale = maximum(sorted_events)-min_t
    T = (sorted_events.-min_t)./scale
    Δt = T[2:end].-T[1:end-1]

    # make sure the events have non-zero waiting times
    if minimum(Δt)<=0.0
        throw(ArgumentError("[KBD ERROR] Events need to be spaced appart by positive times, they cannot occur simultaneously. Preprocess data accordingly."))
    end

    # get the maximum number of states required for optimality
    Kinf = max_state_bound(Δt,s)

    # solve the dynamic program to get estimate for state path
    args_transition = log(n)*γ
    α0 = Float64(n)
    args_data = (α0,s)
    state_sequence, _ = solve_dynamic_program(Δt,kbd_standard_transition_cost,kbd_standard_data_model,Kinf,args_transition,args_data)
   
    # generate the burst hierarch tree structure
    root = hierarchical_burst_structure(sorted_events,state_sequence)
    # since we sorted in the beginning, we shoudl return the sorted events aswell in the result
    result = BurstStructure(sorted_events,state_sequence,root)
    return result
end

end
