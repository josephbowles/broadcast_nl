% This script loads the point (channel and measurements) and the inequality
% which certifies the \eta=0.7355 detection efficiency threshold for the
% maximally entangled two qubit state.
% To make the linear program more robust, instead of doing a 1 shot
% feasibility test, we relax the linear program in the following way.
% Instead of imposing that the probability vector for the hidden variables
% is positive, we impose that it is bigger than lambda, and then maximize
% lambda in the linear program. Thus if the LP returns a positive lamba,
% this means it found a solution to the program and it is feasible, but if
% it gives a negative lamba it means the original problem is infeasible.
% This allows for a smoother transition between feasible and infeasible and
% it is more robust numerically, but it also introduces the uncertainty of
% how negative should lambda be before we consider the problem infeasible.
% This directly affects what efficiency threshold eta we get. Here we use
% the convention of -1e-8. If it is more negative than -1e-8 we consider
% the problem infeasible.

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'aux_funcSec3A3CSec4',filesep);
addpath(newdir2);


ins = [3,2,2];
outs = [2,2,2];
dimscell = num2cell([ins outs]);
nr_parties = length(ins);


%INITIAL_VISIBILITY = 0;  % !! I'm using the convention of (1-p)psi+p*id instead of p*psi+(1-p)*id
%CONST_CHI = 0.155;  % For the partially entangled states
%IDENTITY_PLACEMENT = 'A'; % If using rho_A \otimes Id/2 (IDENTITY_PLACEMENT = 'B') or Id/2 \otimes rho_B (IDENTITY_PLACEMENT = 'B') 
%STATE_SETTINGS = struct('name','partially_entangled','CONST_CHI', CONST_CHI, 'IDENTITY_PLACEMENT', IDENTITY_PLACEMENT);
STATE_SETTINGS = struct('name','werner');

%%

%% Load the point of interest

load('detectionefficiency07355data.mat');

%% Calculate the 'optimizer' objects
fprintf("Starting to calculate the YALMIP optimizer for the linear program.\n");
% tic
% yalmip_optimizer_broadcast_p1p2 = BroadcastLP_optimizer(ins, outs);  % this ones takes as input 2 distribs and gives the feasibility
% toc
tic
yalmip_optimizer_detector_broadcast = BroadcastLPfeasibility_optimizer(ins,outs+1,true);
toc
%tic
%yalmip_optimizer_detector_local = BroadcastLPfeasibility_optimizer(ins,outs+1,false);
%toc
% tic
% yalmip_optimizer_broadcast = BroadcastLPfeasibility_optimizer(ins,outs,true);
% toc
% tic
% yalmip_optimizer_local = BroadcastLPfeasibility_optimizer(ins,outs,false);
% toc

% tic
% optimizer_ch = optimizer_channel(ins, outs, 2);
% toc
% tic
% optimizer_a = optimizer_povm_party(1, ins, outs);
% toc
% tic
% optimizer_b = optimizer_povm_party(2, ins, outs);
% toc
% tic
% optimizer_c = optimizer_povm_party(3, ins, outs);
% toc
% optimizer_objects = {optimizer_ch, optimizer_a, optimizer_b, optimizer_c};

%%

sortedetas = sort(unique(list_of_conv_etas));
best_eta = sortedetas(2); % By eye, the first 1 is 0
best_idx = find(list_of_conv_etas==best_eta);
best_idx = 4658;

channel = list_of_conv_channels{best_idx};
POVMs = list_of_conv_povms{best_idx};

checkThatChannelIsGood(channel,2,4,1e-12);
checkPOVMsAreGood(POVMs,ins,outs,1e-12);
kr = KrausOperators(channel, [2, 4]);

fprintf("Kraus operators of the channel:\n");
disp([kr{:}]);

fprintf("Measurements:\n");
for party=1:nr_parties
    for x=1:ins(party)
        obs_x = POVMs{party}{x}{1}-POVMs{party}{x}{2};
        bloch = BlochComponents(obs_x);
        bloch = num2cell(bloch(2:4));
        [azimuth,elevation,r] = cart2sph(bloch{:});
        azimuth = azimuth*180/pi;
        elevation = elevation*180/pi;
        fprintf("Party: %d, input:%d, obs: (azimuth[º],elevation[º],r)=(%f,%f,%f)\n", party, x, azimuth, elevation, r);
        disp(obs_x);
    end
end

BlochAllObs(POVMs);

%%

p_entangled = ProbMultidimArray(final_state(NoisyState(0, STATE_SETTINGS), channel), POVMs, ins, outs);
p_uniform   = ProbMultidimArray(final_state(NoisyState(1, STATE_SETTINGS), channel), POVMs, ins, outs);

checkThatProbSumsToOne(p_entangled, ins, outs);
checkThatProbSumsToOne(p_uniform, ins, outs);

%[alpha, bellcoeffs, LPstatus, ~] = BroadcastLP(p_entangled, p_uniform);
[alpha, ~, ~, duals] = yalmip_optimizer_broadcast_p1p2(p_entangled, p_uniform); 
duals1 = duals{1};
bellcoeffs = reshape(-duals1(1:numel(p_entangled)), [dimscell{:}]);
fprintf("Visibility of the point with 0.7355 detection efficiency."); 
alpha

bellcoeffs2 = clean(bellcoeffs./min(abs(bellcoeffs(abs(bellcoeffs)>1e-6))),1e-6);
bellcoeffs2corr = vpa(dispBellCoeffsCorrelators(bellcoeffs2, ins, outs)*2,3);

%%


ETA_CONVERGENCE_TOL = 1e-8;
ETA_NEGATIVITY_THRESH = 1e-8;
[new_eta, newbellcoeffs, etaLPstatus] = BroadcastSlowLPdetectorefficiency( ...
    yalmip_optimizer_detector_broadcast,p_entangled,STATE_SETTINGS,channel,POVMs,ins,outs,ETA_CONVERGENCE_TOL,ETA_NEGATIVITY_THRESH);
list_of_conv_etas = [list_of_conv_etas, new_eta];

fprintf("detection efficiency throgh a bisection:\n");
new_eta

fprintf("If lam is positive, the noisy probability distribution has a LHV model. \nIf lam is negative, then it doesn't have a LHV model, thus its outside the corresponding local set.\n");
fprintf("This adds an extra output to the probability distribution and checks memebership in the broadcast local set with 1 extra output.\n");
for eta = flip(0.735:0.0001:0.737)
    noisy_prob3 = giveDetectorEfficiencydistrib(p_entangled, eta * ones(1,nr_parties), ins, outs);
    p_entangled = ProbMultidimArray(final_state(NoisyState(0, STATE_SETTINGS), channel), POVMs, ins, outs);
    %noisy_prob = ProbMultidimArrayDETINEFF(final_state(NoisyState(0, STATE_SETTINGS), channel), POVMs, eta*ones(1,nr_parties), ins, outs);
    %noisy_prob2 = giveDetectorEfficiencydistrib_manual(p_entangled, eta*ones(1,nr_parties), ins, outs);
    %aux_c1=checkThatProbSumsToOne(noisy_prob,ins,outs+1);
    [lam, ~, ~, ~] = yalmip_optimizer_detector_broadcast(noisy_prob3);
    %[~, infeas, ~, ~] = yalmip_optimizer_detector_local(noisy_prob3);
    %[bellcoeffs, infeas] = BroadcastLP_membership(noisy_prob, ins, outs+1, true);
    fprintf("eta=%f lam=%g\n", eta, lam);
end

%%
fprintf("Instead of using 1 extra output, we can also take the Bell inequality \n and calculate its value assuming different determinitic detector responses.\n When optimizing over all possible responses, we recover the same 0.7355 point.\n");
    best_eta1 = 1;
    for lamA = 1:8
        for lamB = 1:4
            for lamC = 1:4
                
                eta1 = 1;
                eta0 = 0;
                mid = (eta1+eta0)/2;
                delta_eta=abs(eta1-eta0);
                
                eta_tol = 1e-8;
                MAX_ITER = 100;
                
                iter_nr = 1;
                while delta_eta>eta_tol && iter_nr<MAX_ITER
                   val1 = eval_eta(bellcoeffs, POVMs, channel, eta1, lamA, lamB, lamC, ins, outs);
                   val0 = eval_eta(bellcoeffs, POVMs, channel, eta0, lamA, lamB, lamC, ins, outs);
                   valmid = eval_eta(bellcoeffs, POVMs, channel, (eta1+eta0)/2, lamA, lamB, lamC, ins, outs);
                   
                   if valmid > 0
                       eta1 = (eta1+eta0)/2;
                   elseif val0 < val1
                       eta0 = (eta1+eta0)/2;
                   end
                    
                   iter_nr = iter_nr + 1;
                end

                val1 = eval_eta(bellcoeffs, POVMs, channel, 1, lamA, lamB, lamC, ins, outs);
                eta = eta1;


                fprintf("detector failure responses: detA %d, detA %d, detA %d, eta: %f\n", lamA, lamB, lamC, eta1);
                
                if eta1 < best_eta1
                    best_eta1 = eta1;
                    fprintf("New best eta1 found: %f\n", best_eta1);
                end
            end
        end
    end
    fprintf("Best detection efficiency threshold after optimizing over all detector failure responses: %f\n", best_eta1);
    %%
    fprintf("The bell inequality in correlator notation:\n");
    bellcorr = dispBellCoeffsCorrelators(bellcoeffs,ins,outs);
    [C,T] = coeffs(bellcorr);
    C=vpa(C,3);
    C(abs(C)<1e-5)=0;
    C=C/min(abs(C(abs(C)>0)));
    bellcorr_clean = vpa(dot(C,T),3);
    disp(bellcorr_clean);
