function [POVMs, channel, finalObj, alpha, FLAG] = give_good_initial_condition_ineq(optimizer_objects, bellcoeffs, localbound, STATE_SETTINGS, VARIOUS_CONSTANTS, ins, outs)
ALPHA_INI_ITER = VARIOUS_CONSTANTS('ALPHA_INI_ITER');
TOL_DIST_TO_POINT = VARIOUS_CONSTANTS('TOL_DIST_TO_POINT');
INITIAL_VISIBILITY = VARIOUS_CONSTANTS('INITIAL_VISIBILITY');
alpha=0;
alpha_iteration=1;
LPstatus = 0;
state_visibility = VARIOUS_CONSTANTS('INITIAL_VISIBILITY');
while(alpha_iteration <= ALPHA_INI_ITER) && (abs(alpha-0)<TOL_DIST_TO_POINT || ...
                                            abs(alpha-1)<TOL_DIST_TO_POINT || ...
                                            alpha < 0 || alpha > 1)
    POVMs = givePprojRANDgeneral(ins);  % Random projective measurements
    channel = {giveChannelRAND(2,4)};  % Random isometry

    %p_entangled = ProbMultidimArray(final_state(NoisyState(INITIAL_VISIBILITY, STATE_SETTINGS), channel), POVMs, ins, outs);
    %p_uniform   = ProbMultidimArray(final_state(NoisyState(1, STATE_SETTINGS), channel), POVMs, ins, outs);
    %[alpha, LPstatus]= visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform);
    %aux1 = bellcoeffs(:)'*p_entangled(:); aux2 = bellcoeffs(:)'*p_uniform(:);
    %fprintf("Vis of ineq given random ini = %f\t ineq_value=%f\t bound_l=%f \n", alpha, aux1, localbound);
    
    %[POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyState(visibility, STATE_SETTINGS), POVMs, channel);
    [POVMs,finalObj,channel] = optimizer_seesaw(optimizer_objects, bellcoeffs, NoisyState(state_visibility, STATE_SETTINGS), POVMs, channel, ins, outs);

    p_entangled = ProbMultidimArray(final_state(NoisyState(INITIAL_VISIBILITY, STATE_SETTINGS), channel), POVMs, ins, outs);
    p_uniform   = ProbMultidimArray(final_state(NoisyState(1, STATE_SETTINGS), channel), POVMs, ins, outs);
    [alpha, LPstatus]= visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform);
    
    aux1 = bellcoeffs(:)'*p_entangled(:); aux2 = bellcoeffs(:)'*p_uniform(:);
    %fprintf("%d/%d Vis of ineq after optimizing = %f\t ineq_value=%f\t bound_l=%f \n", ...
    %        alpha_iteration, ALPHA_INI_ITER, alpha, aux1, localbound);

    if LPstatus == 0 && abs(alpha-0)> TOL_DIST_TO_POINT && abs(alpha-1)> TOL_DIST_TO_POINT &&  abs(finalObj-localbound)>TOL_DIST_TO_POINT
        %aux1 = bellcoeffs(:)'*p_entangled(:); aux2 = bellcoeffs(:)'*p_uniform(:);
        %    fprintf("Vis of ineq given random ini = %f. ineq_value=%f bound_l=%f \n", ...
        %    alpha, aux1, localbound);
        fprintf("   Found good initial condition. Bell value:%f NS2Bound:%f diff=%g ineq_visibility= %f state_vis=%f\n", finalObj, localbound, finalObj-localbound, alpha, state_visibility);
        break;
    else
       fprintf("%d Visibility LP not solver correctly (LPstatus=%d). Trying another initial condition.\n", alpha_iteration, LPstatus);
    end
    alpha_iteration = alpha_iteration + 1;
end
if alpha_iteration > ALPHA_INI_ITER || abs(finalObj-localbound)<TOL_DIST_TO_POINT
   fprintf("Couldn't find a good initial point before hitting the maximum number of iterations. Going to a different inequality.")
   FLAG = false;      
else
   fprintf("Good initial condition. diff=%g\n", finalObj-localbound);
   FLAG = true; 
end

end
