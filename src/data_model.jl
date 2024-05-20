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
