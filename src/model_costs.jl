###
# Data model for the distribution of the waiting times of the events, where parameter s incorporates
# the current activity state of the underlying model.
# A standard model here would be exponential.
# The function give log-likelihood-values
###

function kbd_exponential_waiting_time_model(Δt::Float64,state::Int64,θ::Tuple{Float64,Float64})
    # parameters:
    # Δt: waiting time
    # state: activity level of the automaton
    # θ: tuple of base rate α0 and scaling factor s of the automaton
    α0,s = θ
    α = α0 * s^state
    nL = -state*log(s)+Δt*α
    return nL
end

# set standard to exponential waiting time model
kbd_standard_data_model = kbd_exponential_waiting_time_model

####
# Models for transition cost between automaton states i and j
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
