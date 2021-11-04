%%
% Code to find the visibility of the Chained 33-22 inequality promoted to the
% broadcast scenario with 1 Alice and 2 Bobs. The visibility returned by
% the program uses the opposite convention (vis=0 == maximally entangled
% state and vis=1 == noisy state). After enough iterations it will give as
% best visibility vis=0.61... in the paper (vis=1-0.61..=0.39 in the code).

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'QETLAB-0.9',filesep);
addpath(newdir2);
addpath(genpath(newdir3));

% Fix the seed for debugging.
%rng(1234); % this seems to saturate the local bound
rng('default');

load('bellcoeffsChained3322broadcast.mat');

final_round_alpha = []; % store the final visibility after a round of 
final_round_povm = {};
final_round_channels = {};

% LOOP PARAMETERS
MAX_ITER_OUTER_LOOP = 50;
ALHPA_TOL_DIST_TO_POINT = 1e-3;
DELTA_ALPHA_CONVERGENCE_THRESHOLD = 1e-6;
MAX_ITER_INNER_LOOP = 1000;
ALPHA_STDEV_TOL = 1;
%%%

PSD_TOLERANCE = 1e-8;


% INITIAL VALUES LOOP
latest_alpha_meta = 0;
meta_iteration = 1;

while meta_iteration < MAX_ITER_OUTER_LOOP
    fprintf("\nRound %d of the see-saw.\n", meta_iteration);
    
    [localbound, LPstatus] = ClassicalOptInequality_fromLPBroadcast(bellcoeffs);
    fprintf("Local bound of bell inequality: %f\n", localbound);
    % store outputs from the loop
    inner_list_of_alphas = [];
    inner_list_of_povms = {};
    inner_list_of_channels = {};
    
    % initial values inner loop
    iteration = 1;
    alpha_stdev = ALPHA_STDEV_TOL-1;
    alpha = 0;
    
    while ( iteration<MAX_ITER_INNER_LOOP && ...
            ( abs(alpha-0)<ALHPA_TOL_DIST_TO_POINT || ... % the following two just aim to get a visibility between 0 and 1
            abs(alpha-1)<ALHPA_TOL_DIST_TO_POINT) )

        spv0=-1;
        spv1=0;
        % the point of this lopp is to get AT LEAST get a starting point
        % that isn't worse than just the maximally mixed state
        while spv0<spv1
            %POVMs = givePprojRAND2(ins,outs);
            POVMs = givePprojRANDgeneral(ins);
            %POVMs = givePprojRANDmaxEBI();
            %POVMs = givePprojDET();
            %channel = RandomSuperoperator([2,4]);
            channel = {giveChannelRAND(2,4)};
            %channel = {give_Joe_U()};
            %channel = giveChannelAddsIdentity(2,2,"right");

            spv0 = evaluate_bell_ineq(bellcoeffs, 0, final_state(NoisyWernerState(0), channel), POVMs, ins, outs);
            spv1 = evaluate_bell_ineq(bellcoeffs, 0, final_state(NoisyWernerState(1), channel), POVMs, ins, outs);
            fprintf("s·p(v=0)=%f s·p(v=1)=%f localbound=%f\n", spv0, spv1, localbound);
        end
        
        %p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs);
        p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs, ins, outs);
        p_uniform = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs, ins, outs);
        alpha = visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform);
        fprintf("Starting visibility: %f\n", alpha);

        [POVMs,finalObj,channel] = SeeSawOverAllParties(bellcoeffs, NoisyWernerState(0), POVMs, channel);
        assert(checkPOVMsAreGood(POVMs,ins,outs,PSD_TOLERANCE),'Problem with POVMs');
        assert(checkThatChannelIsGood(channel, 2, 4,PSD_TOLERANCE), 'Problem with the channel');

        p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs, ins, outs);
        p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs, ins, outs);
        alpha = visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform);
        fprintf("Local maximum: %f\n", alpha);

        inner_list_of_alphas = [inner_list_of_alphas, alpha];
        inner_list_of_povms{iteration} = POVMs;
        inner_list_of_channels{iteration} = channel;
        iteration = iteration + 1;
           
    end
    final_round_alpha = [final_round_alpha, inner_list_of_alphas(end)];
    final_round_povm{length(final_round_alpha)} = inner_list_of_povms{end};
    final_round_channels{length(final_round_alpha)} = inner_list_of_channels{end};
    latest_alpha_meta = final_round_alpha(end);
    meta_iteration = meta_iteration + 1;
end

% look at the best overall
best_alpha = max(final_round_alpha);
index = find(final_round_alpha == best_alpha);
best_povm = final_round_povm(index(1));
best_channels = final_round_channels(index(1));

fprintf("Best of %d: %f\n", MAX_ITER_META, best_alpha);

save('scenario3.mat')