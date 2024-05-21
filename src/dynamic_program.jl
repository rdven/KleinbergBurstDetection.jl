
###
# The bound from Kleinbergs paper that reduces the infinite state automaton to a finitie state automaton with 
# the same optimal trajectory estimate.
# the time need to be rescaled to [0,1] domain already
###
function max_state_bound(Δt::Vector{Float64},s::Float64)
    δ = minimum(Δt)
    k = Int(ceil(1-log(δ)/log(s)))
    return k
end

###
# The following function is used to solve the MLE problem as a dynamic program (maximizing log-lieklihood).
# This technically works
###
function solve_dynamic_program(Δt::Vector{Float64},transition_cost,data_model,maxState::Int64,args_transition,args_data)::Tuple{Vector{Int64},Float64}
    n = size(Δt,1)
    # allocate cost matrix, which has length (n,maxState+1) we allow any starting state with no penalization
    COST = zeros(Float64,(n,maxState+1))
    PTR = zeros(Int64,(n,maxState+1))
    # inital state has data cost and transition cost starting from 0
    maxInd = maxState+1
    for sind ∈ 1:maxInd
        COST[1,sind] = data_model(Δt[1],sind-1,args_data)
        COST[1,sind] += transition_cost(0,sind-1,args_transition)
    end
    # populate the transition cost matrix
    TC = zeros(Float64,(maxState+1,maxState+1))
    for s1 in 1:maxInd
        for s2 in 1:maxInd
            TC[s1,s2] = transition_cost(s1-1,s2-1,args_transition)
        end
    end
    # now solve the dynamic programming by populating the cost matrix
    for t in 2:n
        for sind ∈ 1:maxInd
            # data cost
            COST[t,sind] = data_model(Δt[t],sind-1,args_data)
            # identify best previous+transiation cost option
            best_prev = 1
            best_cost = COST[t-1,1] + TC[1,sind]
            for sprev ∈ 2:maxInd
                ctemp = COST[t-1,sprev] + TC[sprev,sind]
                if ctemp<best_cost
                    best_prev = sprev
                    best_cost = ctemp
                end
            end
            COST[t,sind] += best_cost
            PTR[t,sind] = best_prev
        end
    end
    # now find the best end state
    best_end = 1
    best_final_cost = COST[n,1]
    for sfinal ∈ 2:maxInd
        if COST[n,sfinal]<best_final_cost
            best_end = sfinal
            best_final_cost = COST[n,sfinal]
        end
    end
    opt_trajectory = zeros(Int64,n)
    opt_trajectory[end] = best_end
    for tminus ∈ 1:(n-1)
        current_t = n-tminus
        opt_trajectory[current_t] = PTR[current_t+1,opt_trajectory[1+current_t]]
    end
    # reduce to account for the index
    opt_trajectory = opt_trajectory .- 1
    return opt_trajectory,best_final_cost
end
