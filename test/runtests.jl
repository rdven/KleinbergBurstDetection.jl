using KleinbergBurstDetection
using Test

function get_test_data_set_1()
    return [0.0,1.0,2.0,2.1,2.2,2.3,2.4,2.5,3.0,4.0,5.0,6.0,6.1,6.11,6.12,6.13,6.14,6.15,6.2,6.3,6.4,6.41,6.42,6.43,6.44,6.45,6.5,7.0,8.0,9.0,10.0]
end

function get_test_data_set_2()
    return [352.4,353.9229398148148,353.9416782407407,353.9782638888889,353.9784027777778,354.4848842592593,354.52557870370373,354.5262847222222,354.54912037037036,354.5921412037037,354.59465277777775,355.0350462962963,355.0695949074074,355.72541666666666,355.9331712962963,356.40060185185183,356.40150462962964,356.78372685185184,360.0184375,360.0298148148148,360.0416319444444,360.0416435185185,360.0422222222222,360.55922453703704,360.689525462963,360.71747685185187,360.71760416666666,360.75789351851853,360.77510416666667,360.9102314814815,360.91261574074076,361.86599537037034,361.89664351851854,362.00553240740743,362.11894675925925,363.7627662037037,365.18943287037035,365.1894675925926,365.64131944444443,365.72444444444443,365.72447916666664,365.8759259259259,366.97664351851853,367.4893055555556,367.48935185185184,367.5799421296296,367.69460648148146,367.6946412037037,367.9602662037037,368.1533912037037,368.1534027777778,368.1534375,368.5750462962963,368.5750810185185,368.7336921296296,368.7337268518518,368.7959490740741,368.8039699074074,370.1917013888889]
end

function test_upper_bound_states(data,s=2.0)
    T = data.-minimum(data)
    T ./= maximum(T)
    Δt = T[2:end]-T[1:end-1]
    return KleinbergBurstDetection.max_state_bound(Δt,s)
end

@testset "KleinbergBurstDetection.jl" begin

    # tests for the state bound for dynamic programming
    @test test_upper_bound_states(get_test_data_set_1(),2.0) == 11
    @test test_upper_bound_states(get_test_data_set2_(),2.0) == 25

    #####
    # test result on dataset 2
    ####
    result_data_2 = detect_bursts(get_test_data_set_2(),2.0,1.0)

    ### we expect 3 different main bursts
    @test size(result_data_2.hierarchy.children,1) == 3

    
end
