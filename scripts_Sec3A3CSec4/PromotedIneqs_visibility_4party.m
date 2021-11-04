%%
% Code to find the visibilities of inequalities promoted to the broadcast
% scenario with 2 alices and 2 bobs. The promoted inequalities are CHSH,
% the Chained inequality with 3 inputs and the Elegant Bell inequality.
% The visibilities are the same as in the 1 alice 2 bobs broadcast
% scenario, namely 0.5777, 0.61 and 0.68. Note that the code uses the
% opposite convention and it returns 1 minus the above visibilities.

%%

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'aux_funcSec3A3CSec4',filesep);
addpath(newdir2);


%% CHOOSE HERE THE INEQUALITY

%load('bellcoeffs_chained_4party_broadcast.mat');
%load('bellcoeffs_ebi_4party_broadcast.mat');
load('bellcoeffs_chsh_4party_broadcast.mat');

%%
bellcoeffs = bellcoeffs;

% ins = size(bellcoeffs,1:4);
% outs = size(bellcoeffs,5:8);
% 
% corr2bell = fromCorrToBellMat(ins,outs);
% 
% bellcoeffs = reshape(corr2bell * corr_mat(2:end)',[ins outs]);
% 

%%

list_of_vis = [];

localbound = double(2*localbound);

for Iteration=1:100
    
    fprintf("\n NEW ROUND OF INITIAL CONDITIONS \n");

POVMs = givePprojRANDgeneral(ins);
%disp(BlochComponents(POVMs{1}{1}{1}))

channels = cell(2,1);
channels{1} = {giveChannelRAND(2,4)};
channels{2} = {giveChannelRAND(2,4)};
% channels{1} = RandomSuperoperator([2,4]);
% channels{2} = RandomSuperoperator([2,4]);

state0 = final_state2(NoisyWernerState(0), channels{1}, channels{2});
state1 = final_state2(NoisyWernerState(1), channels{1}, channels{2});

p_e = CalcProbArray(final_state2(NoisyWernerState(0), channels{1}, channels{2}), POVMs, ins,outs);
p_u = CalcProbArray(final_state2(NoisyWernerState(1), channels{1}, channels{2}), POVMs, ins,outs);

%disp(norm(p_e(:)-p_u(:)))
%disp(norm(final_state2(NoisyWernerState(0), channels{1}, channels{2})-final_state2(NoisyWernerState(1), channels{1}, channels{2})))

fprintf("%f %f %f\n", p_e(:)'*bellcoeffs(:), p_u(:)'*bellcoeffs(:), localbound);

CONV_TOL = 1e-6;
MAX_ITER = 50;
obj_val = -1e6;
delta_obj = 1e6;



iter = 1;
while delta_obj > CONV_TOL || iter > MAX_ITER
   
    
    [newChoiMap,newObjective,problemStatus] = SeeSawOverChannel2(NoisyWernerState(0), bellcoeffs, POVMs, channels, 'A', ins, outs);
    channels{1} = newChoiMap;
    fprintf("After optimizing channel over Alice, new objective=%f\n", newObjective); 
    delta_obj = abs(newObjective-obj_val);
    obj_val = newObjective;
%     if delta_obj < CONV_TOL
%         break; 
%     end
    p_e = CalcProbArray(final_state2(NoisyWernerState(0), channels{1}, channels{2}), POVMs, ins,outs);
%     fprintf("Obje=%f\n", sum(p_e.*bellcoeffs,'all'));
    state0 = final_state2(NoisyWernerState(0), channels{1}, channels{2});
    
        [newChoiMap,newObjective,problemStatus] = SeeSawOverChannel2(NoisyWernerState(0), bellcoeffs, POVMs, channels, 'B', ins, outs);
    channels{2} = newChoiMap;
    fprintf("After optimizing channel over Bob, new objective=%f\n", newObjective); 
    delta_obj = abs(newObjective-obj_val);
    obj_val = newObjective;
%     if delta_obj < CONV_TOL
%         break; 
%     end
    p_e = CalcProbArray(final_state2(NoisyWernerState(0), channels{1}, channels{2}), POVMs, ins,outs);
%     fprintf("Obje=%f\n", sum(p_e.*bellcoeffs,'all'));
    state0 = final_state2(NoisyWernerState(0), channels{1}, channels{2});
    

    
    
        [POVMs,newObjective,problemStatus] = SeeSawOverASingleParty2(3, state0, bellcoeffs, POVMs, ins, outs);
    fprintf("After optimizing POVMs over p=%d, new objective=%f\n", 3, newObjective); 
    delta_obj = abs(newObjective-obj_val);
    obj_val = newObjective;
%     if delta_obj < CONV_TOL
%         break; 
%     end
%     p_e = CalcProbArray(final_state2(NoisyWernerState(0), channels{1}, channels{2}), POVMs, ins,outs);
%     fprintf("Obje=%f\n", sum(p_e.*bellcoeffs,'all'));
    
    [POVMs,newObjective,problemStatus] = SeeSawOverASingleParty2(4, state0, bellcoeffs, POVMs, ins, outs);
    fprintf("After optimizing POVMs over p=%d, new objective=%f\n", 4, newObjective); 
    delta_obj = abs(newObjective-obj_val);
    obj_val = newObjective;
%     if delta_obj < CONV_TOL
%         break; 
%     end
    
%     p_e = CalcProbArray(final_state2(NoisyWernerState(0), channels{1}, channels{2}), POVMs, ins,outs);
%     fprintf("Obje=%f\n", sum(p_e.*bellcoeffs,'all'));


    [POVMs,newObjective,problemStatus] = SeeSawOverASingleParty2(1, state0, bellcoeffs, POVMs, ins, outs);
    fprintf("After optimizing POVMs over p=%d, new objective=%f\n", 1, newObjective); 
    delta_obj = abs(newObjective-obj_val);
    obj_val = newObjective;
%     if delta_obj < CONV_TOL
%         break; 
%     end
%     p_e = CalcProbArray(final_state2(NoisyWernerState(0), channels{1}, channels{2}), POVMs, ins,outs);
%     fprintf("Obje=%f\n", sum(p_e.*bellcoeffs,'all'));
    
    [POVMs,newObjective,problemStatus] = SeeSawOverASingleParty2(2, state0, bellcoeffs, POVMs, ins, outs);
    fprintf("After optimizing POVMs over p=%d, new objective=%f\n", 2, newObjective); 
        delta_obj = abs(newObjective-obj_val);
    obj_val = newObjective;
%     if delta_obj < CONV_TOL
%         break; 
%     end
%     p_e = CalcProbArray(final_state2(NoisyWernerState(0), channels{1}, channels{2}), POVMs, ins,outs);
%     fprintf("Obje=%f\n", sum(p_e.*bellcoeffs,'all'));
    

    
%     if obj_val > localbound
%         fprintf("Nice! We got a point outside.\n");
%     end

    p_e = CalcProbArray(final_state2(NoisyWernerState(0), channels{1}, channels{2}), POVMs, ins,outs);
    p_u = CalcProbArray(final_state2(NoisyWernerState(1), channels{1}, channels{2}), POVMs, ins,outs);
    be = p_e(:)'*bellcoeffs(:);
    bu = p_u(:)'*bellcoeffs(:);
    visibility = (localbound-bu)/(be-bu);
%     fprintf("visibility = %f (%f %f)\n", , be, bu);
%     

    iter = iter + 1;
end
if iter > MAX_ITER
   warning('Hit maximum number of iterations!'); 
end
fprintf("Iter: %d Converged to: %f (localbound: %f) visibility: %f (%f %f)\n", Iteration, newObjective, localbound, visibility, be, bu); 
list_of_vis = [list_of_vis, visibility];

fprintf("Best visibilities:\n");
sortedvis = sort(list_of_vis);
if size(sortedvis)>10
   disp(sortedvis(1:10));
else
   disp(sortedvis); 
end


end

%%


