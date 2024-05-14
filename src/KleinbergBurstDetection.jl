module KleinberBurstDetection

###
# The following function is used to solve the MLE problem as a dynamic program (maximizing log-lieklihood).
# This technically works
###
function solveDynamicProgram(transition_cost,data_cost,maxState::Int64)::Tuple{Vector{Int64},Float64}
    # TODO
end

###
# The bound from kleinbergs paper that reduces the infinite state automaton to a finitie state automaton with 
# the same optimal trajectory estimate.
###
function maxStateBound(event_times::Vector{Float64},s::Float64)
    T = maximum(event_times)-minimum(event_times)
    Δ = event_times[2:end]-event_times[1:end-1]
    δ = minimum(Δ)
    k = Int(ceil(1+log(T/δ)/log(s)))
    return k 
end

function kbd_standard_transition_cost()

end

function kbd_standard_data_cost()

end

function burstDetection(event_times::Vector{Float64},s=2.0)
    # check values
    if s <= 1.0
        throw(ArgumentError("[KBD ERROR] the scaling factor s needs to be larger than 1.0."))
    end
    if size(event_times,1) < 3
        throw(ArgumentError("[KBD ERROR] there need to be at least 3 events in the time series."))
    end
    # get the maximum number of states required
    Kinf = maxStateBound(event_times,s)
    # define the costs


    # solve the dynamic program to get estimate for state path
    X = solveDynamicProgram(t_cost,d_cost,Kinf)
end

end