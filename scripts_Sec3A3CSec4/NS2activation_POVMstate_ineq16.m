% This script optimizes over a channel and the party measurements to find
% an example of NS genuinte network nonlocality activation. It uses
% inequality 16 from arxiv 1112.2626 but on line 40 on can optimize over
% any inequality desired. On line 39 one can set over which angles chi to
% look at. After enough time, it should converge to p=0.182642, which means
% that for p<=0.182642 there is an inequality certifying the state is NS2
% nonlocal, but for p>=0.155505 the original bipartite state has an LHV
% model for POVMs.

%clear all

diary aux_chidep.txt

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'aux_funcSec3A3CSec4',filesep);
addpath(newdir);
addpath(newdir2);


%% Scenario settings
load('bellcoeffs_arxiv1112_2626.mat'); % loads 'bellcoeffs_cell','local_upper_bounds','ins','outs'
load('table3_arXiv1112_2626.mat'); % loads 'table3arXiv11122626'

ins = [2,2,2];
outs = [2,2,2];
nrparties = length(ins);

%% Calculate the 'optimizer' objects
optimizer_ch = optimizer_channel(ins, outs, 2);
optimizer_a = optimizer_povm_party(1, ins, outs);
optimizer_b = optimizer_povm_party(2, ins, outs);
optimizer_c = optimizer_povm_party(3, ins, outs);
optimizer_objects = {optimizer_ch, optimizer_a, optimizer_b, optimizer_c};

%% 
NR_OF_INEQS = size(bellcoeffs_cell,2);
results_per_ineq = cell(1, 3);
results_within_loop = cell(1,3);

tic
chivalues = [0.05];
chiresults = cell(length(chivalues),NR_OF_INEQS,9);
for ineq_nr=[16]%,10,11,12,13,14,15,16,18,21,26,27,31,32,35,37,41,42,44,51,56,59,73,89,92]%save%NR_OF_INEQS
for chi_idx = 1:length(chivalues)
    % Set seed randomly so that we can reproduce the statistics within one big loop
    %rng('shuffle','twister');
    %fprintf("rng info:\n\n");
    %disp(rng);
    
    bellcoeffs = bellcoeffs_cell{ineq_nr};
    localboundNS2 = local_upper_bounds(ineq_nr);
    quantumbound = table3arXiv11122626(ineq_nr, 5);
    
    fprintf("\n\tInequality number = %d\n", ineq_nr);
    
    % Loop parameters
    ALPHA_INI_ITER = 50; % How many initial conditions to try for the initial point
    MAX_ITER_BIG_LOOP = 15; % 4How many times to try the whole process
    MAX_ITER_VIS_OPT_LOOP = 100; % given good initial condition, loop for trying to optimize the visibility
    
    INITIAL_VISIBILITY = 0;  % !! I'm using the convention of (1-p)psi+p*id instead of p*psi+(1-p)*id
    CONST_CHI = chivalues(chi_idx);%0.05;  % For the partially entangled states
    IDENTITY_PLACEMENT = 'A'; % If using rho_A \otimes Id/2 (IDENTITY_PLACEMENT = 'B') or Id/2 \otimes rho_B (IDENTITY_PLACEMENT = 'A') 
    %STATE_SETTINGS = struct('name','partially_entangled','CONST_CHI', CONST_CHI, 'IDENTITY_PLACEMENT', IDENTITY_PLACEMENT);
    STATE_SETTINGS = struct('name','pent_povm','CONST_CHI', CONST_CHI, 'IDENTITY_PLACEMENT', IDENTITY_PLACEMENT);
    %STATE_SETTINGS = struct('name','werner');
    
    % Solve ineq (23) from https://arxiv.org/pdf/1510.06721.pdf
    COS2 = (cos(2*CONST_CHI))^2;
    roots23 = roots([COS2,-2*COS2,0,2,-1]);
    roots23 = roots23(abs(imag(roots23))<1e-8); % discard imaginary roots;
    roots23 = roots23(abs(roots23)>=0);
    roots23 = roots23(abs(roots23)<=1); % discard p outside [0,1]
    roots23 = max(roots23); % just in case there is more than one, but there shouldnt be
    p = roots23 - 0.1; % check that ineq (23) is satisfied in the intervat (0,roots23)
    assert(dot([COS2,-2*COS2,0,2,-1],[p^4,p^3,p^2,p^1,1]) <= 0, "Something bad happened");
    % for all p >
    UNSTEERABILITY_THRESHOLD = 1 - roots23; % we're using opposite convention for noise
    
    TOL_DIST_TO_POINT = 1e-4; % Util for excluding "bad" points such 
                                    % as visibility=0 or visibility=1 (if 
                                    % its zero, then most likely its an infeasibility)
    ALPHA_CONVERGENCE_THRESHOLD = 1e-6;
    VISIBILITY_CONVERGENCE_THRESHOLD = 1e-6;
    
    DELTA_STATE_VIS_PROP = 0.95;

    VARIOUS_CONSTANTS = containers.Map();
    VARIOUS_CONSTANTS('ALPHA_INI_ITER') = ALPHA_INI_ITER;
    VARIOUS_CONSTANTS('INITIAL_VISIBILITY') = INITIAL_VISIBILITY;
    VARIOUS_CONSTANTS('TOL_DIST_TO_POINT') = TOL_DIST_TO_POINT;
    VARIOUS_CONSTANTS('VISIBILITY_CONVERGENCE_THRESHOLD') = VISIBILITY_CONVERGENCE_THRESHOLD;
    VARIOUS_CONSTANTS('DELTA_STATE_VIS_PROP') = DELTA_STATE_VIS_PROP;
    VARIOUS_CONSTANTS('MAX_ITER_VIS_OPT_LOOP') = MAX_ITER_VIS_OPT_LOOP;
    
    
    % Calculate the best per inequality
    best_visibility = 0;
    best_visibility_iter = 1;
    best_POVMS = givePprojRANDgeneral(ins);
    best_channel = {giveChannelRAND(2,4)};
    best_visisibility_ineq_nr = 1;
    results_within_loop = cell(1);
    state_visibility = 0;
    
    last_good_channel = 0;
    last_good_POVMs = 0;
    
    iteration_big_loop = 1;
    while iteration_big_loop <= MAX_ITER_BIG_LOOP 
        [POVMs, channel, finalObj, alpha, FLAG] = give_good_initial_condition_ineq(optimizer_objects, bellcoeffs, localboundNS2, STATE_SETTINGS, VARIOUS_CONSTANTS, ins, outs);
        if FLAG == false
            ("Couldn't find a good initial condition! Moving to a different inequality.");
            state_visibility = 0;
            break;
        else
            [POVMs, channel, state_visibility] = optimize_state_visibility(optimizer_objects, bellcoeffs, localboundNS2, POVMs, channel, STATE_SETTINGS, VARIOUS_CONSTANTS, ins, outs);
        end
        fprintf("%d Final state visibility       : \t%f\n", iteration_big_loop, state_visibility);

            
        if state_visibility > best_visibility
            best_visibility = state_visibility;
            best_visibility_iter = iteration_big_loop;
            best_POVMS = last_good_POVMs;
            best_channel = last_good_channel;
            best_visisibility_ineq_nr = ineq_nr;
            
            if strcmp(STATE_SETTINGS.name, 'pent_povm')
                if STATE_SETTINGS.IDENTITY_PLACEMENT == 'B'
                    fprintf("\n\tNew best visibility found! p=%f (UNS A->B, Steerable B->A for p>=%f, i.e., LHV for p>=%f)\n", state_visibility, UNSTEERABILITY_THRESHOLD, UNSTEERABILITY_THRESHOLD);
                elseif STATE_SETTINGS.IDENTITY_PLACEMENT == 'A'
                    fprintf("\n\tNew best visibility found! p=%f (UNS B->A, Steerable A->B for p>=%f, i.e., LHV for p>=%f)\n", state_visibility, UNSTEERABILITY_THRESHOLD, UNSTEERABILITY_THRESHOLD);
                end
            elseif strcmp(STATE_SETTINGS.name, 'werner')
                fprintf("\n\tNew best visibility found! p=%f\n", state_visibility);
            end
                
%             fprintf("Channel:\n");
%             disp(best_channel);
%             fprintf("Channel spectrum=\n");
%             disp(eig(best_channel).');
%             fprintf("Measurements:\n");
%             for party=1:nrparties
%                 for x=1:ins(party)
%                     obs_x = POVMs{party}{x}{1}-POVMs{party}{x}{2};
%                     bloch = BlochComponents(obs_x);
%                     bloch = num2cell(bloch(2:4));
%                     [azimuth,elevation,r] = cart2sph(bloch{:});
%                     azimuth = azimuth*180/pi;
%                     elevation = elevation*180/pi;
%                     fprintf("Party: %d, input:%d, obs: (azimuth[º],elevation[º],r)=(%f,%f,%f)\n", party, x, azimuth, elevation, r);
%                 end
%             end
           
        end
        results_within_loop{iteration_big_loop, 1} = state_visibility;
        
        iteration_big_loop = iteration_big_loop + 1;
    end
    all_vis_for_ineq = [results_within_loop{:,1}];
    promedio = mean(all_vis_for_ineq);
    desviacion = sqrt(var(all_vis_for_ineq));
    
    results_per_ineq{ineq_nr, 1} = best_visibility;
    results_per_ineq{ineq_nr, 2} = best_channel;
    results_per_ineq{ineq_nr, 3} = best_POVMS;
    results_per_ineq{ineq_nr, 4} = promedio;
    results_per_ineq{ineq_nr, 5} = desviacion; 
    results_per_ineq{ineq_nr, 6} = [results_within_loop{:,1}]; 
    
    all_vis = [results_per_ineq{:,1}];
    wheremax = all_vis == max(all_vis);
    auxints = 1:length(all_vis);
    position_max = auxints(wheremax); % this might be degenerate i just take the first
    value_max = all_vis(wheremax);
    
    
    
    %fprintf("\n\t Current best is state_visibility=%f for inequality %d. Average: %g +- %g (1 sigma).\n\n", value_max(1), position_max(1), promedio, desviacion);
    

    % Save to file
    fprintf("Saving current workspace to file.\n");
    ScenarioFilename = 'scenario';
    filename = strcat(ScenarioFilename,'aux_chi_dep','.mat');
    save(filename);
    
fprintf("\nCHI LOOP chi bestfinalstatevisibility pthresh ineqnr = %f %f %f %d\n\n\n", CONST_CHI, value_max(1), UNSTEERABILITY_THRESHOLD, ineq_nr);
chiresults{chi_idx, ineq_nr, 1} = CONST_CHI;
chiresults{chi_idx, ineq_nr, 2} = value_max(1);
chiresults{chi_idx, ineq_nr, 3} = UNSTEERABILITY_THRESHOLD;
chiresults{chi_idx, ineq_nr, 4} = best_visibility;
chiresults{chi_idx, ineq_nr, 5} = best_channel;
chiresults{chi_idx, ineq_nr, 6} = best_POVMS;
chiresults{chi_idx, ineq_nr, 7} = bellcoeffs;
chiresults{chi_idx, ineq_nr, 8} = localboundNS2;
chiresults{chi_idx, ineq_nr, 9} = quantumbound;
   

end

end
toc
%save('aux_chi_dep_after8oct.mat');

%%
% hold on
% 
% %yline(1/sqrt(2),'-','','LineWidth',1,'DisplayName','1/√2 line');
% ineq_nrs = [6,7,8,9];
% ineq_idx = ineq_nrs(1);
% 
% zdata = [chiresults{:,ineq_idx,3}];
% plot([chiresults{:,ineq_idx,1}],1-zdata,'Color','r','DisplayName','p_{threshold}','LineWidth',1.8)
% 
% white_c = [1 1 1];
% black_c = [0 0 0];
% 
% idx = 1;
% for i=ineq_nrs
%     xdata = [chiresults{:,i,1}];
%     ydata = [chiresults{:,i,2}];
%     idx = idx + 1;
%     plot(xdata,1-ydata,'DisplayName',num2str(i),'LineWidth',1.5);%,'Color',0.5*white_c*(1-idx/length(ineq_nrs))/10+(0.8*ones(1,3))*(1-idx/length(ineq_nrs)),'LineWidth',1);%,'DisplayName','p_{critical}','LineWidth',1.8)
% end
% 
% best_vis = zeros(size(chiresults,1),1);
% for i=1:size(chiresults,1)
%     best_vis(i) = max([chiresults{i,:,2}]);
% end
% %plot([chiresults{:,ineq_idx,1}],1-best_vis,'Color',[0 0 0],'DisplayName','p_{threshold}','LineWidth',1.5)
% 
% 
% title('Visibility upper bounds for all $NS_2$ inequalities', 'Interpreter', 'latex');
% xlim([0 1]);
% ylim([0 1]);
% xlabel('$\chi$','Interpreter', 'latex')
% ylabel('$p$', 'Interpreter', 'latex')
% 
% xticks([0 0.3 0.6 pi/4 0.9 1.2 pi/2, 1.8]);
% xticklabels({'0','0.3','0.6','π/4','0.9','1.2','π/2','1.8'});
% legend('Location','southeast')
% 
% %annotation('textbox', [0.279+0.175 0.215 0.215 0.], 'interpreter','latex', 'EdgeColor','w','String','$\rho_{POVM}=\frac{1}{2} \rho(p,\chi) + \frac{1}{2} \rho_A(p,\chi) \otimes \vert 0 \rangle\langle 0\vert$','FitBoxToText','on')
% 
% %grid on
% 
% 
% hold off


%%
%plot([results_per_ineq{:,1}])

function [POVMs,newObj,channel] = optimizer_seesaw(optimizer_objects, bellcoeffs, state_AB, POVMs, channel, ins, outs)

optimizer_ch = optimizer_objects{1};
optimizer_a = optimizer_objects{2};
optimizer_b = optimizer_objects{3};
optimizer_c = optimizer_objects{4};

%% Loop parameters
MAX_NR_ITERATIONS = 500;
CONVERGENCE_TOL = 1e-6;
%%

deltaObj = 1e6;
oldObj = -1e6;
newObj = 0;

iteration = 1;
while deltaObj > CONVERGENCE_TOL && iteration <= MAX_NR_ITERATIONS 
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
    
    newObj = output{1};
    deltaObj = abs(newObj-oldObj);

    iteration = iteration + 1;
end
if iteration > MAX_NR_ITERATIONS
   warning("Hit maximum number of iterations (=%d) without conv. with tol=%g. Returning last value as converged.", MAX_NR_ITERATIONS, CONVERGENCE_TOL); 
end

end

