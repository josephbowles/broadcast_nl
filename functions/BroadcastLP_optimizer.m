function optim_object = BroadcastLP_optimizer(nr_inputs_per_party, nr_outputs_per_party)
    % IMPORTANT alpha is defined as (1-alpha)*p1 + (alpha)*p2
    % usually p1 is the entangled, and p2 the uniform, so this
    % program minimizes alpha using this convention
    % if p1 and p2 are both local then all alpha are good
    % if p1 and p2 are both nonlocal the program is infeasible
    % ideally p1 should be nonlocal

    % IMPORTANT I use the convention where the probabilities are called as
    % p(x,y,z,a,b,c) instead of p(a,b,c,x,y,z)
%     dims_p = size(p1);
%     
%     assert(mod(length(dims_p),2) == 0, "The probability array should have equal number of inputs and outputs.");
%     assert(all(size(p1) == size(p2)), "The two probability arrays should have equal dimenisons.");
%     if ( norm(p1(:) - p2(:)) < 1e-2)
%         warning( "The probabilities should not be almost equal. norm(p1-p2)=%g", norm(p1(:) - p2(:)));
%         alpha=0;
%         bellcoeffs=0;
%         dual_alpha=0;
%         LPstatus = 1;
%     end
    
    nrparties = length(nr_inputs_per_party);
    
    %nr_inputs_per_party = dims_p(1:nrparties);
    %nr_outputs_per_party = dims_p(nrparties+1:end);

    party_for_det_points = 1; % this is party 'A' % TODO Make this a function input
    nrinputsofA  = nr_inputs_per_party(party_for_det_points);
    nroutputsofA = nr_outputs_per_party(party_for_det_points);
    
    nr_det_points = nroutputsofA^nrinputsofA;
  
    alpha = sdpvar(1); % alpha will be the visibility

    tempdims = [nr_det_points, nr_inputs_per_party(2:end), nr_outputs_per_party(2:end)];
    q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
    tempdims_cell = num2cell(tempdims);
    qarray = sdpvar(prod(tempdims(:)),1);
    q = cell(tempdims_cell{:});
    for idx = 1:size(q_coords,1)
        coords = num2cell(q_coords(idx,:));
        q{coords{:}} = qarray(idx);
        allins = fullfact(nr_inputs_per_party(2:end));
        allouts = fullfact(nr_outputs_per_party(2:end));
        normalisationconstraint = [];
        for x_idx = 1:size(allins,1)
            input_slice = num2cell(allins(x_idx,:));
            suma = 0;
            for a_idx = 1:size(allouts,1)
                output_slice = num2cell(allouts(a_idx,:));
                for lam = 1:nr_det_points
                    suma = suma + q{lam,input_slice{:},output_slice{:}};
                end
            end
            normalisationconstraint = [normalisationconstraint, suma == 1];
        end
    end
    
    
    %q = sdpvar(nr_det_points, inputs_cell{2:end}, outputs_cell{2:end});
    
    visibility_constraints = [alpha >= 0];

%    positivityconstraints = [];
% 	auxsize=size(q);
%     for i=1:prod(auxsize(:))
%         positivityconstraints = [positivityconstraints, qarray(i) >= 0];
%     end
    positivityconstraints = [qarray(:) >= 0];
    
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

    dimscell = num2cell([nr_inputs_per_party nr_outputs_per_party]);
    p1 = sdpvar(dimscell{:});
    p2 = sdpvar(dimscell{:});
    noisy_prob = reshape((1-alpha) * p1(:) + (alpha) * p2(:), [dimscell{:}]);
    
    probability_constraints = [];
    productstructure = [nr_inputs_per_party, nr_outputs_per_party];
    cartesianproduct_forprobconstraints = ind2subv(productstructure, 1:prod(productstructure(:)));
    for idx = 1:size(cartesianproduct_forprobconstraints,1)
        coords_cell = num2cell(cartesianproduct_forprobconstraints(idx,:));
        summ = 0;
        for lam = 1:nr_det_points
            % for 3 parties, same as det_strategy(lam,x,a) * q(lam,y,z,b,c);
            summ = summ + det_strategy(lam, coords_cell{1}, coords_cell{nrparties+1}) ... 
                * q{lam, coords_cell{2:nrparties}, coords_cell{nrparties+2:end}};
        end
        probability_constraints = [probability_constraints, summ == noisy_prob(coords_cell{:})];
    end
    
    %% Solving the SDP
    
    objective = alpha;
            
    constraints = [probability_constraints,...
        normalisationconstraint, ...
        positivityconstraints, ...
        visibility_constraints, ...
        nonsignalling_constraintsB, ...
        nonsignalling_constraintsC, ...
        nonsignalling_constraintsBC];
    %optsol = optimize(constraints, objective, ...
    %    sdpsettings('solver','mosek','verbose',0,'showprogress',0, 'debug', 0, 'warning',0));
        
    params_in = {p1,p2};
    sols_out = {alpha};
    opts = sdpsettings('solver','mosek','verbose',0);
    
    %optsol = optimize(constraints, objective, ...
    %    sdpsettings('solver','mosek','verbose',0,'showprogress',0, 'debug', 0, 'warning',0));
    optim_object = optimizer(constraints, objective, opts, params_in, sols_out);
    
%                       'mosek.MSK_DPAR_INTPNT_TOL_INFEAS', 1e-6,...
%                       'mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS',1e-6));
%     optsol = optimize(constraints, objective, ...
%         sdpsettings('solver','gurobi','verbose',1,'showprogress',1, 'debug', 1, 'warning',1));    
%   
%    LPstatus = optsol.problem;
%     if LPstatus ~= 0
%         disp(optsol)
%         error('Check why the problem is not successfully solved.');
%     end
    
%     dimcell = num2cell([nr_inputs_per_party,nr_outputs_per_party]);
%     bellcoeffs = zeros(dimcell{:});           
%     index = 1;
%     for idx = 1:size(cartesianproduct_forprobconstraints,1)
%         % Here I need to make sure that say, index = 6, corresponds to the
%         % same (x,y,z,a,b,c) tuple that it did when defining the
%         % probability constraints. I make sure of that by using the same
%         % array for looping, 'cartesianproduct_forprobconstraints'
%         coords_cell = num2cell(cartesianproduct_forprobconstraints(idx,:));
%         bellcoeffs(coords_cell{:}) = -value(dual(probability_constraints(index)));
%         index = index + 1;
%     end
%     final_alpha = value(alpha);
%     dual_alpha = 0;
    %fprintf("lam4 · (p1-p2) = %f\n", sum(bellcoeffs .* (p1-p2),'all'));
end
