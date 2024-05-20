####
# Models for transiation cost between automaton states i and j
####

# the transition costs proposed by Kleinberg
function kbd_standard_transition_cost(i::Int64,j::Int64,θ::Float64)::Float64
    # the extra parameter here is just a control variable how much penalized switching states is
    # Kleinberg sets θ = ln(n)*γ where n is the number of inter event time intervals
    if j <= i
        return 0.0
    else
        return θ*(j-i)
    end
end
