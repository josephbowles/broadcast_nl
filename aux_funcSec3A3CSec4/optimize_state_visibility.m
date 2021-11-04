function [POVMs, channel, visibility] = optimize_state_visibility(optimizer_objects, bellcoeffs, localboundNS2, POVMs, channel, STATE_SETTINGS, VARIOUS_CONSTANTS, ins, outs)

MAX_ITER_VIS_OPT_LOOP = VARIOUS_CONSTANTS('MAX_ITER_VIS_OPT_LOOP');
VISIBILITY_CONVERGENCE_THRESHOLD = VARIOUS_CONSTANTS('VISIBILITY_CONVERGENCE_THRESHOLD');
TOL_DIST_TO_POINT = VARIOUS_CONSTANTS('TOL_DIST_TO_POINT');
DELTA_STATE_VIS_PROP = VARIOUS_CONSTANTS('DELTA_STATE_VIS_PROP');

state_visibility = 0;

p_entangled = ProbMultidimArray(final_state(NoisyState(state_visibility, STATE_SETTINGS), channel), POVMs, ins, outs);
p_uniform   = ProbMultidimArray(final_state(NoisyState(1,                STATE_SETTINGS), channel), POVMs, ins, outs);
[ineq_visibility, LPstatus]= visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);
if abs(ineq_visibility-0)<TOL_DIST_TO_POINT || ineq_visibility < 0 || ineq_visibility > 1
   warning("Not a good initial point!"); 
end
%fprintf("Initial inequality visibility: \t%f\n", ineq_visibility);

iter_vis_opt_loop = 1;
%old_visibility = -1;
AbsDeltaStateVis = 1e6;
while iter_vis_opt_loop <= MAX_ITER_VIS_OPT_LOOP && AbsDeltaStateVis > VISIBILITY_CONVERGENCE_THRESHOLD
    
    old_visibility = state_visibility;
    state_visibility = state_visibility + ineq_visibility*DELTA_STATE_VIS_PROP - state_visibility*ineq_visibility*DELTA_STATE_VIS_PROP; % remember I use the opposite convention for visibility
    AbsDeltaStateVis = abs(state_visibility-old_visibility);
    %fprintf("State visibility: %f %d\n", state_visibility, iter_vis_opt_loop);
    [POVMs,~,channel] = optimizer_seesaw(optimizer_objects, bellcoeffs, NoisyState(state_visibility, STATE_SETTINGS), POVMs, channel, ins, outs);
    %channel = cleanChannel(channel, 2, 4);
    %checkPOVMsAreGood(POVMs,ins,outs);
    %checkThatChannelIsGood(channel, 2, 4);

    p_entangled = ProbMultidimArray(final_state(NoisyState(state_visibility, STATE_SETTINGS), channel), POVMs, ins, outs);
    p_uniform   = ProbMultidimArray(final_state(NoisyState(1,                STATE_SETTINGS), channel), POVMs, ins, outs);
    [ineq_visibility, LPstatus]= visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);


    if abs(ineq_visibility-1)<TOL_DIST_TO_POINT || ineq_visibility > 1 || ineq_visibility < 0 || LPstatus ~= 0
        warning("Getting ridiculous visibilities (vis=%f LPstatus=%d). Trying different initial conditions.\n");
        [POVMs, channel, finalObj, alpha, FLAG] = give_good_initial_condition_ineq(optimizer_objects, ...
                 bellcoeffs, localboundNS2, STATE_SETTINGS, VARIOUS_CONSTANTS, ins, outs);
        % WARNING: ASSUMING THAT FLAG == true, else you shouldn't get
        % inside this function to optimize the visibility 
    end

    if abs(ineq_visibility-0) > TOL_DIST_TO_POINT
        % This is because if the best ineq_visibility is 0,
        % then probably the output channel and POVMs from the
        % optimization are thrash, so we should keep the last
        % non trivial ones.
        last_good_channel = channel;
        last_good_POVMs = POVMs;
    end
    
    iter_vis_opt_loop = iter_vis_opt_loop + 1;
end

%fprintf("Final state visibility       : \t%f\n", state_visibility);
if iter_vis_opt_loop > MAX_ITER_VIS_OPT_LOOP
    warning("Couldn't converge to any good visibility within MAX ITER.");
end
    visibility = state_visibility;

end

