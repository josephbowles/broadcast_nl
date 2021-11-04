function [eta, bellcoeffs, LPstatus] = BroadcastSlowLPdetectorefficiency(LP_feas_optimizer,ini_prob,STATE_SETTINGS,channel,POVMs,ins,outs,ETA_CONV_TOLERANCE,TOL_NEGATIVITY_BISECTION)

% for this code p1 should be outside the local set and p2 inside

%TOL_NEGATIVITY_BISECTION = 1e-7;

dims_p = [ins(:)', outs(:)'];
aux_dims_detector_eff = num2cell([ins, outs + 1]);
aux_cell = num2cell(dims_p);
eta = 0;
eta_m = 0;
bellcoeffs = zeros(aux_cell{:});

nr_parties = length(ins);

%p1 = giveDetectorEfficiencydistrib_manual(ini_prob, 1.0 * ones(1,nr_parties), ins, outs);
%p2 = giveDetectorEfficiencydistrib_manuaL1l(ini_prob, 0.0 * ones(1,nr_parties), ins, outs);
p1 = ProbMultidimArrayDETINEFF(final_state(NoisyState(0, STATE_SETTINGS), channel), POVMs, 1*ones(1,nr_parties), ins, outs);
p2 = ProbMultidimArrayDETINEFF(final_state(NoisyState(0, STATE_SETTINGS), channel), POVMs, 0*ones(1,nr_parties), ins, outs);


[lamval1, infeas1] = LP_feas_optimizer(p1);%SlowLP(p1, nr_inputs_per_party, nr_outputs_per_party);
[lamval2, infeas2] = LP_feas_optimizer(p2);%SlowLP(p2, nr_inputs_per_party, nr_outputs_per_party);
%fprintf("%g %g\n", lamval1, lamval2);

if infeas1 == true || infeas2 == true
   eta = 0;
    aux_cell = num2cell(dims_p);
    bellcoeffs = zeros(aux_cell{:});
    LPstatus = false;
    warning("One of the LP's is infeasible!");
    return;
end

% if L1 == true == 1, then it means the LP is infeasible
% if L1 == false == 0, then it means the LP is feasible
if lamval1 >= -TOL_NEGATIVITY_BISECTION && lamval2 >= -TOL_NEGATIVITY_BISECTION  % both are feasible, therefore inside the sets
    eta = 0;
    aux_cell = num2cell(dims_p);
    bellcoeffs = zeros(aux_cell{:});
    LPstatus = false;
    warning("Both points are inside the local set!");
    return;
end
if lamval1 <= -TOL_NEGATIVITY_BISECTION && lamval2 <= -TOL_NEGATIVITY_BISECTION  % both are infeasible, therefore outside
    eta = 0;
    aux_cell = num2cell(dims_p);
    bellcoeffs = zeros(aux_cell{:});
    LPstatus = false;
    warning("Both points are outside the local set!");
    return;
end
if lamval1 >= -TOL_NEGATIVITY_BISECTION && lamval2 <= -TOL_NEGATIVITY_BISECTION  % it's the wrong order, p1 should be infeasible
    eta = 0;
    aux_cell = num2cell(dims_p);
    bellcoeffs = zeros(aux_cell{:});
    LPstatus = false;
    warning("Better to switch around p1 and p2 (so p1 is outside the broadcast-local set)");
    return;
end
if lamval1 <= -TOL_NEGATIVITY_BISECTION && lamval2 >= -TOL_NEGATIVITY_BISECTION
    % this algorithm is assuming that the function is monotone, thus
    % there is a single critical eta where the the function goes from
    % feasible to infeasible
    
    eta1 = 1;  % infeasible point, lamval < 0
    eta2 = 0;  % feasible point, lamval > 0
    eta_m = 1/2;  % unknown
    
    %infeas_eta1 = true;  % P(eta1) is infeasible
    %infeas_eta2 = false;  % P(eta2) is feasible
    noisy_prob = ProbMultidimArrayDETINEFF(final_state(NoisyState(0, STATE_SETTINGS), channel), POVMs, eta_m*ones(1,nr_parties), ins, outs);
    [lamval, infeas_i] =  LP_feas_optimizer(noisy_prob); %SlowLP(noisy_prob, nr_inputs_per_party, nr_outputs_per_party + 1);

    precision = 1e6;
    while precision>ETA_CONV_TOLERANCE && infeas_i == false
        precision = abs((eta2-eta1)/2);
        if lamval < -TOL_NEGATIVITY_BISECTION %infeas_eta_m == true  % eta_m is infeasible
            eta1 = eta_m;
        else
            eta2 = eta_m;
        end
        eta_m = (eta1+eta2)/2; 
        
        noisy_prob =  ProbMultidimArrayDETINEFF(final_state(NoisyState(0, STATE_SETTINGS), channel), ...
            POVMs, eta_m*ones(1,nr_parties), ins, outs);
        [lamval, infeas_i] =  LP_feas_optimizer(noisy_prob); %SlowLP(noisy_prob, nr_inputs_per_party, nr_outputs_per_party + 1);
        %fprintf("LP bisection progress: am=%g lm=%g\n", eta_m, infeas_eta_m);
       
    end

%     if infeas_eta_m == true
%         [~, infeas_eta_m, ~, bellcoeffs_temp] = LP_feas_optimizer(noisy_prob);
%         bellcoeffs_temp1 = bellcoeffs_temp{1};
%         bellcoeffs = reshape(bellcoeffs_temp1(1:prod([aux_dims_detector_eff{:}])),[aux_dims_detector_eff{:}]);  % TODO CHECK CORRECTNESS
%     else
%         noisy_prob = giveDetectorEfficiencydistrib(ini_prob, (eta_m+precision)*ones(1,nr_parties), ins, outs);
%         [~, infeas_eta_m, ~, bellcoeffs_temp] = LP_feas_optimizer(noisy_prob);
%         bellcoeffs_temp1 = bellcoeffs_temp{1};
%         bellcoeffs = reshape(bellcoeffs_temp1(1:prod([aux_dims_detector_eff{:}])),[aux_dims_detector_eff{:}]);  % TODO CHECK CORRECTNESS
%     end
    fprintf("lamval=%g infeas_i=%f\n", lamval, infeas_i);
    LPstatus = true;
end

eta = eta_m;
% bellcoeffs is already calculated
end

function [LPstatus, bellcoeffs] = SlowLP(noisy_prob, nr_inputs_per_party, nr_outputs_per_party, flag_broadcast_set)
    dims_p = [nr_inputs_per_party(:)', nr_outputs_per_party(:)'];
    
    nrparties = length(nr_inputs_per_party);
    
    %nr_inputs_per_party = dims_p(1:nrparties);
    %nr_outputs_per_party = dims_p(nrparties+1:end);

    party_for_det_points = 1; % this is party 'A' % TODO Make this a function input
    nrinputsofA  = nr_inputs_per_party(party_for_det_points);
    nroutputsofA = nr_outputs_per_party(party_for_det_points);
    
    nr_det_points = nroutputsofA^nrinputsofA;
  

    tempdims = [nr_det_points, nr_inputs_per_party(2:end), nr_outputs_per_party(2:end)];
    q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
    tempdims_cell = num2cell(tempdims);
    qarray = sdpvar(prod(tempdims(:)),1);
    q = cell(tempdims_cell{:});
    for idx = 1:size(q_coords,1)
        coords = num2cell(q_coords(idx,:));
        q{coords{:}} = qarray(idx);
    end

    positivityconstraints = [];
	auxsize=size(q);
    for i=1:prod(auxsize(:))
        positivityconstraints = [positivityconstraints, qarray(i) >= 0];
    end
    
    
    %% Non signalling constraints
    
    % Note that we are imposing that q(bd|yzw lam)=q(bd|y w lam) (no z
    % dependence) but we don't yet marginalize over d. This is because if
    % the previous is true, then clearly \sum_d q(bd|yw lam) = q(b| y w
    % lam) is not dependent on z. The condition "q(bd|yzw lam) indep of z" is
    % stronger than "q(b|yzw lam) indep of z"
    
    nonsignalling_constraintsBD = [];
    % non signaling for bob:
    for lam = 1:nr_det_points
        coordstructure = [nr_outputs_per_party(2), nr_inputs_per_party(2), nr_inputs_per_party(4), nr_outputs_per_party(4)];
        all_b_y_and_w_d = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_b_y_and_w_d,1)
            b = all_b_y_and_w_d(slice,1);
            y = all_b_y_and_w_d(slice,2);
            w = all_b_y_and_w_d(slice,3);
            d = all_b_y_and_w_d(slice,4);
            % choose z = 1 for summ1
            summ1 = 0;
            for c = 1:nr_outputs_per_party(3)
               summ1 = summ1 + q{lam,y,1,w,b,c,d};
            end
            % the marginal for z != 1 should be equal to z=1
            for z = 2:nr_inputs_per_party(3)
                summ2 = 0;
                for c = 1:nr_outputs_per_party(3)
                   summ2 = summ2 + q{lam,y,z,w,b,c,d};
                end
                nonsignalling_constraintsBD = [nonsignalling_constraintsBD, summ1 == summ2];
            end
        end
    end

    %non signaling for charlie:
    nonsignalling_constraintsCD = [];
    for lam = 1:nr_det_points
        coordstructure = [nr_outputs_per_party(3), nr_inputs_per_party(3), nr_inputs_per_party(4), nr_outputs_per_party(4)];
        all_c_z_and_w_d = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_c_z_and_w_d,1)
            c = all_c_z_and_w_d(slice,1);
            z = all_c_z_and_w_d(slice,2);
            w = all_c_z_and_w_d(slice,3);
            d = all_c_z_and_w_d(slice,4);
            % choose y = 1
            summ1 = 0;
            for b = 1:nr_outputs_per_party(2)
               summ1 = summ1 + q{lam,1,z,w,b,c,d};
            end
            % the marginal for z' != 1 should be equal to z=1
            for y=2:nr_inputs_per_party(2)
                summ2 = 0;
                for b = 1:nr_outputs_per_party(2)
                   summ2 = summ2 + q{lam,y,z,w,b,c,d};
                end
                nonsignalling_constraintsCD = [nonsignalling_constraintsCD, summ1 == summ2];
            end
        end
    end
    
    %non signaling between the instrument and bob and charlie
    % nonsignaling for instrument
    % nonsignaling D not influenced by B1 and B2 (or B and C)
    nonsignalling_constraintsD = [];
     for lam = 1:nr_det_points
        coordstructure = [nr_outputs_per_party(4), nr_inputs_per_party(4)];
        all_d_and_w = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for idx = 1:size(all_d_and_w,1)
            d = all_d_and_w(idx,1);
            w = all_d_and_w(idx,2);
            
            summ1 = 0;
            for b = 1:nr_outputs_per_party(2)
                for c = 1:nr_outputs_per_party(3)
                    summ1 = summ1 + q{lam,1,1,w,b,c,d};
                end
            end
            
            inputstructure = [nr_inputs_per_party(2), nr_inputs_per_party(3)];
            YZinputsCartesianProduct = ind2subv(inputstructure, 1:prod(inputstructure(:)));
            for index = 2:size(YZinputsCartesianProduct,1)
                y2 = YZinputsCartesianProduct(index,1);
                z2 = YZinputsCartesianProduct(index,2);
                summ2 = 0;
                for b = 1:nr_outputs_per_party(2)
                    for c = 1:nr_outputs_per_party(3)
                        summ2 = summ2 + q{lam,y2,z2,w,b,c,d};
                    end
                end
                nonsignalling_constraintsD = [nonsignalling_constraintsD, ...
                                        summ1 == summ2];
            end
        end
     end
     
    % partial normalization nonsignaling (summing over bcd should not depend
    % on the inputs)   
    nonsignalling_constraintsBCD = [];
    for lam = 1:nr_det_points          
        inputstructure = [nr_inputs_per_party(2), nr_inputs_per_party(3), nr_inputs_per_party(4)];
        all_y_and_z_and_w = ind2subv(inputstructure, 1:prod(inputstructure(:)));
        
        slice = 1;
        y1 = all_y_and_z_and_w(slice,1);
        z1 = all_y_and_z_and_w(slice,2);
        w1 = all_y_and_z_and_w(slice,3);
        summ1 = 0;
        for b = 1:nr_outputs_per_party(2)
            for c = 1:nr_outputs_per_party(3)
                for d = 1:nr_outputs_per_party(4)
                    summ1 = summ1 + q{lam,y1,z1,w1,b,c,d};
                end
            end
        end
        
        for slice = 2:size(all_y_and_z_and_w,1)
            y2 = all_y_and_z_and_w(slice,1);
            z2 = all_y_and_z_and_w(slice,2);
            w2 = all_y_and_z_and_w(slice,3);
            summ2 = 0;
            for b = 1:nr_outputs_per_party(2)
                for c = 1:nr_outputs_per_party(3)
                    for d = 1:nr_outputs_per_party(4)
                        summ2 = summ2 + q{lam,y2,z2,w2,b,c,d};
                    end
                end
            end
            nonsignalling_constraintsBCD = [nonsignalling_constraintsBCD, ...
                                                summ1 == summ2];
        end
    end
    
    %% Probability constraints
    det_strategy = givedetstratA(nr_outputs_per_party(1),nr_inputs_per_party(1));
    
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
    
    objective = 0;
            
    constraints = [positivityconstraints, ...
                    probability_constraints, ...
                    nonsignalling_constraintsBD, ...
                    nonsignalling_constraintsCD, ...
                    nonsignalling_constraintsD, ...
                    nonsignalling_constraintsBCD];
    optsol = optimize(constraints, objective, ...
        sdpsettings('solver','mosek','verbose',0));

    LPstatus = optsol.problem;
    
    dimcell = num2cell([nr_inputs_per_party,nr_outputs_per_party]);
    bellcoeffs = zeros(dimcell{:});           
    
    if LPstatus == 0
        for index = 1:size(cartesianproduct_forprobconstraints,1)
            % Here I need to make sure that say, index = 6, corresponds to the
            % same (x,y,z,a,b,c) tuple that it did when defining the
            % probability constraints. I make sure of that by using the same
            % array for looping, 'cartesianproduct_forprobconstraints'
            coords_cell = num2cell(cartesianproduct_forprobconstraints(index,:));
            bellcoeffs(coords_cell{:}) = -value(dual(probability_constraints(index)));
        end
    end

end

