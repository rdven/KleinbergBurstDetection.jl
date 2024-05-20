module KleinbergBurstDetection

include("data_model.jl")
include("transition_costs.jl")
include("dynamic_program.jl")

###
# The bound from Kleinbergs paper that reduces the infinite state automaton to a finitie state automaton with 
# the same optimal trajectory estimate.
# the time needs to be rescaled to [0,1] domain already
###
function maxStateBound(Δt::Vector{Float64},s::Float64)
    δ = minimum(Δt)
    k = Int(ceil(1-log(δ)/log(s)))
    return k
end


function detect_bursts(event_times::Vector{Float64},s=2.0,γ=1.0)
    # check values
    if s <= 1.0
        throw(ArgumentError("[KBD ERROR] the scaling factor s needs to be larger than 1.0."))
    end
    if size(event_times,1) < 3
        throw(ArgumentError("[KBD ERROR] there need to be at least 3 events in the time series."))
    end
    # rescale event times
    n = size(event_times,1)
    min_t = minimum(event_times)
    scale = maximum(event_times)-min_t
    T = (event_times.-min_t)./scale
    Δt = T[2:end].-T[1:end-1]
    # make sure the events have non-zero waiting times
    if minimum(Δt)==0.0
        throw(ArgumentError("[KBD ERROR] Events need to be spaced appart by positive times, they cannot occur simultaneously. Preprocess data accordingly."))
    end
    # get the maximum number of states required
    Kinf = maxStateBound(Δt,s)
    # solve the dynamic program to get estimate for state path
    args_transition = log(n)*γ
    α0 = Float64(n)
    args_data = (α0,s)
    state_path, _ = solveDynamicProgram(Δt,kbd_standard_transition_cost,kbd_standard_data_model,Kinf,args_transition,args_data)
    return state_path
end

end
