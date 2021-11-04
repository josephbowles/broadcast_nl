function [POVMs,newObj,channel] = optimizer_seesaw(optimizer_objects, bellcoeffs, state_AB, POVMs, channel, ins, outs)

optimizer_ch = optimizer_objects{1};
optimizer_a = optimizer_objects{2};
optimizer_b = optimizer_objects{3};
optimizer_c = optimizer_objects{4};
% yalmip_optimizer_broadcast_p1p2 = optimizer_objects{5};

%% Loop parameters
MAX_NR_ITERATIONS = 500;
CONVERGENCE_TOL = 1e-6;
%%

deltaObj = 1e6;
oldObj = -1e6;
newObj = 0;
LPstatus = 0;
iteration = 1;
while deltaObj > CONVERGENCE_TOL && iteration <= MAX_NR_ITERATIONS  && LPstatus == 0
    oldObj = newObj;
    %% Channel
    belloperator = give_Bell_operator(bellcoeffs, POVMs, ins, outs);
    state = state_AB;
    ia_state = give_ia_state(state);

    %output = optimizer_ch([{ia_state}, {belloperator}]);
    output = optimizer_ch({ia_state,belloperator});
    channel = output{1};
    
%     newObj = output{2};

    %% Do the POVMs
    output_state = final_state(state, channel);
    ia_output_state = give_ia_state(output_state);

    %% Alice
    partial_products_for_a = give_partial_products(POVMs, bellcoeffs, 1, ins, outs);
    output = optimizer_a([{ia_output_state}, partial_products_for_a(:)']);
%     oldObj = newObj;
%     newObj = output{1};
%     deltaObj = abs(newObj-oldObj);
    povm_alice = reshape({output{2:end}},[ins(1),outs(1)]);
    for x=1:ins(1)
        for a=1:outs(1)
            POVMs{1}{x}{a} = povm_alice{x,a};  % Update the POVMs
        end
    end

    %% Bob
    partial_products_for_b = give_partial_products(POVMs, bellcoeffs, 2, ins, outs);
    output = optimizer_b([{ia_output_state}, partial_products_for_b(:)']);
%     oldObj = newObj;
%     newObj = output{1};
%     deltaObj = abs(newObj-oldObj);
    povm_bob = reshape({output{2:end}},[ins(2),outs(2)]);
    for y=1:ins(2)
        for b=1:outs(2)
            POVMs{2}{y}{b} = povm_bob{y,b};  % Update the POVMs
        end
    end

    %% Charlie
    partial_products_for_c = give_partial_products(POVMs, bellcoeffs, 3, ins, outs);
    output = optimizer_c([{ia_output_state}, partial_products_for_c(:)']);
    
    povm_charlie = reshape({output{2:end}},[ins(3),outs(3)]);
    for z=1:ins(3)
        for c=1:outs(3)
            POVMs{3}{z}{c} = povm_charlie{z,c};  % Update the POVMs
        end
    end
    
%     p_entangled = ProbMultidimArray(final_state(NoisyState(0, STATE_SETTINGS), channel), POVMs, ins, outs);
%     p_uniform   = ProbMultidimArray(final_state(NoisyState(1, STATE_SETTINGS), channel), POVMs, ins, outs);
%     [alpha, LPstatus, ~, duals] = yalmip_optimizer_broadcast_p1p2(p_entangled, p_uniform); 
%     duals1 = duals{1};
%     old_bell = bellcoeffs;
%     bellcoeffs = reshape(-duals1(1:numel(p_entangled)), [ins outs]);
%     fprintf("iter alpha LPstatus = %d %f %d\n",iteration, alpha, LPstatus);
    
    if LPstatus == 0
    %newObj = alpha;
    newObj = output{1};
    deltaObj = abs(newObj-oldObj);
    end
    
    iteration = iteration + 1;
end
if iteration > MAX_NR_ITERATIONS
   warning("Hit maximum number of iterations (=%d) without conv. with tol=%g. Returning last value as converged.", MAX_NR_ITERATIONS, CONVERGENCE_TOL); 
end

end
