function optim_object = BroadcastLPfeasibility_optimizer(nr_inputs_per_party, nr_outputs_per_party, flag_broadcast_set)
    nrparties = length(nr_inputs_per_party);

    party_for_det_points = 1; % this is party 'A' % TODO Make this a function input
    nrinputsofA  = nr_inputs_per_party(party_for_det_points);
    nroutputsofA = nr_outputs_per_party(party_for_det_points);
    
    nr_det_points = nroutputsofA^nrinputsofA;
    
    aux_dims = num2cell([nr_inputs_per_party, nr_outputs_per_party]);
    prob = sdpvar(aux_dims{:});
    
    if flag_broadcast_set
        tempdims = [nr_det_points, nr_inputs_per_party(2:end), nr_outputs_per_party(2:end)];
        q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
        tempdims_cell = num2cell(tempdims);
        qarray = sdpvar(prod(tempdims(:)),1);
        q = cell(tempdims_cell{:});
        for idx = 1:size(q_coords,1)
            coords = num2cell(q_coords(idx,:));
            q{coords{:}} = qarray(idx);
        end
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
    else
        nrdetpoints_allparties = num2cell(nr_outputs_per_party.^nr_inputs_per_party);
        %nrhiddenvars = prod(nrdetpoints_allparties(:));
        qarray = sdpvar(nrdetpoints_allparties{:},'full');
        normalisationconstraint = [sum(qarray(:)) == 1];
    end
    
    
    %q = sdpvar(nr_det_points, inputs_cell{2:end}, outputs_cell{2:end});
    
    lampos = sdpvar(1);
    positivityconstraints = [qarray(:) >= lampos];   
    
if flag_broadcast_set
    % If this flag is set, then we check membership to the set defined
    % by p(abc|xyz)=\sum_lam q_lam D(a|x\lam) p^NS(bc|yz \lam), \sum_lam
    % q_lam = 1 and q_lam >= 0
    % where p^NS(bc|yz \lam) is no-signaling and we further do the trick of
    % absorbing q_lam into p^NS(bc|yz \lam).
    % If the flag is set to false, then this is simply the 3-partite local
    % set:
    % p(abc|xyz)=\sum_lam q_lam D(a|x\lam) D(b|y\lam) D(c|z\lam)
    % \sum q_lam = 1 and q_lam >= 0
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
else
    % if we test membership in the local set, we leave these empty
    nonsignalling_constraintsB = [];
    nonsignalling_constraintsC = [];
    nonsignalling_constraintsBC = [];
end
    
    %% Probability constraints

if flag_broadcast_set
    
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
        probability_constraints = [probability_constraints, summ == prob(coords_cell{:})];
    end
else
    %nrdetpoints_allparties = nr_outputs_per_party.^nr_inputs_per_party;
    all_lam_combs = fullfact([nrdetpoints_allparties{:}]);
    all_in_out_combs = fullfact([nr_inputs_per_party, nr_outputs_per_party]);
   
    detA = givedetstratA(nr_outputs_per_party(1),nr_inputs_per_party(1));
    detB = givedetstratA(nr_outputs_per_party(2),nr_inputs_per_party(2));
    detC = givedetstratA(nr_outputs_per_party(3),nr_inputs_per_party(3));

    probability_constraints = [];
       for xa_slice = 1:size(all_in_out_combs,1)
            coord = all_in_out_combs(xa_slice,:);
            summ = 0;
            q_idx = 1;
            for lam_slice = 1:size(all_lam_combs,1)
                lamA = all_lam_combs(lam_slice,1);
                lamB = all_lam_combs(lam_slice,2);
                lamC = all_lam_combs(lam_slice,3);
                detprod = detA(lamA,coord(1),coord(nrparties+1))* ...
                          detB(lamB,coord(2),coord(nrparties+2))* ...
                          detC(lamC,coord(3),coord(nrparties+3));
                      
                %fprintf("%f %f\n", q_idx, size(all_lam_combs,1))
                if detprod ~= 0
                    summ = summ + qarray(q_idx)*detprod;
                end    
                q_idx = q_idx + 1;
            end
            coordcell = num2cell(coord);
            probability_constraints = [probability_constraints, summ == prob(coordcell{:})];
       
       end
    
end
    %% Solving the SDP
    
    objective = -lampos;
            
    constraints = [ probability_constraints, ...
                    positivityconstraints, ...
                    normalisationconstraint, ...
                    nonsignalling_constraintsB, ...
                    nonsignalling_constraintsC, ...
                    nonsignalling_constraintsBC];
    
    params_in = {prob};
    sols_out = {lampos};
    opts = sdpsettings('solver','mosek','verbose',0);
    
    %optsol = optimize(constraints, objective, ...
    %    sdpsettings('solver','mosek','verbose',0,'showprogress',0, 'debug', 0, 'warning',0));
    optim_object = optimizer(constraints, objective, opts, params_in, sols_out);

%    LPstatus = optsol.problem;
    
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
end
