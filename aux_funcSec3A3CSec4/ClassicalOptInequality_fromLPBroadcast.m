function [local_bound, LPstatus] = ClassicalOptInequality_fromLPBroadcast(bellcoeffs) 
    % IMPORTANT I use the convention where the probabilities are called as
    % p(x,y,z,a,b,c) instead of p(a,b,c,x,y,z)
    dims_p = size(bellcoeffs);
    
    assert(mod(length(dims_p),2) == 0, "The bellcoeffs array should have equal number of inputs and outputs.");

    nrparties = length(dims_p)/2;
    
    nr_inputs_per_party = dims_p(1:nrparties);
    nr_outputs_per_party = dims_p(nrparties+1:end);

    party_for_det_points = 1; % this is party 'A' % TODO Make this a function input
    nrinputsofA  = nr_inputs_per_party(party_for_det_points);
    nroutputsofA = nr_outputs_per_party(party_for_det_points);
    
    nr_det_points = nroutputsofA^nrinputsofA;
  
    alpha = sdpvar(1); % alpha will be the visibility

    tempdims = [nr_det_points, nr_inputs_per_party(2:end), nr_outputs_per_party(2:end)];

    q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
    tempdims_cell = num2cell(tempdims);
    qarray = sdpvar(prod(tempdims(:)));
    q = cell(tempdims_cell{:});
    for idx = 1:size(q_coords,1)
        coords = num2cell(q_coords(idx,:));
        q{coords{:}} = qarray(idx);
    end
    
    visibility_constraints = [alpha >= 0];
    %visibility_constraints = [visibility_constraints,alpha <= 1];
    
    positivityconstraints = [];
	auxsizeq=size(q);
    for i=1:prod(auxsizeq(:))
        positivityconstraints = [positivityconstraints, qarray(i) >= 0];
    end
    
    
    %% Non signalling constraints
    
    nonsignalling_constraintsB = [];
    % non signaling for bob:
    for lam = 1:nr_det_points
        coordstructure = [nr_outputs_per_party(2), nr_inputs_per_party(2)];
        all_b_and_y = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_b_and_y,1)
            b = all_b_and_y(slice,1);
            y = all_b_and_y(slice,2);
            % choose z = 1
            summ1 = 0;
            for c = 1:nr_outputs_per_party(3)
                summ1 = summ1 + q{lam,y,1,b,c};
            end
            % the marginal for z != 1 should be equal to z=1
            for z=2:nr_inputs_per_party(3)
                summ2 = 0;
                for c = 1:nr_outputs_per_party(3)
                    summ2 = summ2 + q{lam,y,z,b,c};
                end
                nonsignalling_constraintsB = [nonsignalling_constraintsB, summ1 == summ2];
            end
        end
    end

    %non signaling for charlie:
    nonsignalling_constraintsC = [];
    for lam = 1:nr_det_points
        coordstructure = [nr_outputs_per_party(3), nr_inputs_per_party(3)];
        all_b_and_y = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_b_and_y,1)
            c = all_b_and_y(slice,1);
            z = all_b_and_y(slice,2);
            % choose y = 1
            summ1 = 0;
            for b = 1:nr_outputs_per_party(2)
                summ1 = summ1 + q{lam,1,z,b,c};
            end
            % the marginal for z' != 1 should be equal to z=1
            for y=2:nr_inputs_per_party(2)
                summ2 = 0;
                for b = 1:nr_outputs_per_party(2)
                    summ2 = summ2 + q{lam,y,z,b,c};
                end
                nonsignalling_constraintsC = [nonsignalling_constraintsC, summ1 == summ2];
            end
        end
    end

    % partial normalization nonsignaling (summing over bc should not depend
    % on the inputs)   
    nonsignalling_constraintsBC = [];
    for lam = 1:nr_det_points          
        inputstructure = [nr_inputs_per_party(2), nr_inputs_per_party(3)];
        all_y_and_z = ind2subv(inputstructure, 1:prod(inputstructure(:)));
        
        slice = 1;
        y1 = all_y_and_z(slice,1);
        z1 = all_y_and_z(slice,2);
        summ1 = 0;
        for b = 1:nr_outputs_per_party(2)
            for c = 1:nr_outputs_per_party(3)
                summ1 = summ1 + q{lam,y1,z1,b,c};
            end
        end
        
        for slice = 2:size(all_y_and_z,1)
            y2 = all_y_and_z(slice,1);
            z2 = all_y_and_z(slice,2);
            summ2 = 0;
            for b = 1:nr_outputs_per_party(2)
                for c = 1:nr_outputs_per_party(3)
                    summ2 = summ2 + q{lam,y2,z2,b,c};
                end
            end
            nonsignalling_constraintsBC = [nonsignalling_constraintsBC, summ1 == summ2];
        end
    end
    
    %% Probability constraints
    det_strategy = givedetstratA(nr_outputs_per_party(1),nr_inputs_per_party(1));

    objective = 0;
    normalization_constraints = [];
    %productstructure = [nr_inputs_per_party, nr_outputs_per_party];
    settings_combinations = ind2subv(nr_inputs_per_party, 1:prod(nr_inputs_per_party(:)));
    outputs_combinations = ind2subv(nr_outputs_per_party, 1:prod(nr_outputs_per_party(:)));
    %cartesianproduct_forprobconstraints = ind2subv(productstructure, 1:prod(productstructure(:)));
    for slice_settings = 1:size(settings_combinations,1)
        xyz_cell = num2cell(settings_combinations(slice_settings,:));
        normalization_summ = 0;
        for slice_outputs = 1:size(outputs_combinations, 1)
            abc_cell = num2cell(outputs_combinations(slice_outputs,:));
            summ = 0;
            for lam = 1:nr_det_points
                % for 3 parties, same as det_strategy(lam,x,a) * q(lam,y,z,b,c);
                summ = summ + det_strategy(lam, xyz_cell{1}, abc_cell{1}) ... 
                            * q{lam, xyz_cell{2:end}, abc_cell{2:end}};
            end
            objective = objective + bellcoeffs(xyz_cell{:}, abc_cell{:}) * summ;
            normalization_summ = normalization_summ + summ;
        end
        normalization_constraints = [normalization_constraints, normalization_summ == 1];
    end
    %% Solving the SDP
            
    constraints = [positivityconstraints, ...
        visibility_constraints, ...
        normalization_constraints, ...
        nonsignalling_constraintsB, ...
        nonsignalling_constraintsC, ...
        nonsignalling_constraintsBC];
    
    optsol = optimize(constraints, -objective, ...
        sdpsettings('solver','mosek','verbose',0,'showprogress',0, 'debug', 0, 'warning',0));
    
%                       'mosek.MSK_DPAR_INTPNT_TOL_INFEAS', 1e-6,...
%                       'mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS',1e-6));
%     optsol = optimize(constraints, objective, ...
%         sdpsettings('solver','gurobi','verbose',1,'showprogress',1, 'debug', 1, 'warning',1));    
%   
    LPstatus = optsol.problem;
    if LPstatus ~= 0
        disp(optsol)
        error('Check why the problem is not successfully solved.');
    end    
    local_bound = value(objective);
end

