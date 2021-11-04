function funcoutput = ClassicalOptInequalityNoBroadcast(bellcoeffs) 
    assert(mod(length(size(bellcoeffs)),2)==0,"There should be as many inputs as outputs.");
    dims = size(bellcoeffs);
    nrparties = length(dims)/2;
    ins = dims(1:nrparties);
    outs = dims(nrparties+1:end);
    
    nr_det_points = [];
    for party=1:nrparties
       nr_det_points = [nr_det_points, outs(party)^ins(party)]; 
    end

    q = sdpvar(prod(nr_det_points),1);
                         
    positivityconstraints = [q >= 0];

    cartproductOUT = ind2subv(outs, 1:prod(outs(:)));
    cartproductIN  = ind2subv(ins,  1:prod(ins(:)));
    
    dets = {};
    for party=1:nrparties
       dets{party} = givedetstratA(outs(party),ins(party));
    end
    
    objective = 0;
    norm_constraints = [];
   for i2 = 1:size(cartproductIN,1)
       settings = num2cell(cartproductIN(i2,:));
       normalization_summ = 0;
       for i1 = 1:size(cartproductOUT,1)
          outputs = num2cell(cartproductOUT(i1,:));
          summ = 0;
          for lam = 1:prod(nr_det_points)
              [lam_comb] = ind2subv(nr_det_points,lam);
              auxprod = 1;
              for party=1:nrparties
                 auxarray = dets{party};
                 auxprod = auxprod * auxarray(lam_comb(party),settings{party},outputs{party});
              end
              summ = summ + q(lam) * auxprod;
          end
          objective = objective + bellcoeffs(settings{:},outputs{:})*summ;
          normalization_summ = normalization_summ + summ;
       end
       norm_constraints = [norm_constraints, normalization_summ == 1];
   end
    
    constraints = [positivityconstraints, norm_constraints];

    optsol = optimize(constraints, -objective, ...
        sdpsettings('solver','mosek', ...
    'verbose',0,'dualize',0, ...
    'showprogress',0,...
    'savesolverinput',0,'savesolveroutput',0,'debug',0));

% uncomment if you have probablity constraints and you want the bell
% inequality
%     nrduals = length(probability_constraints);
%     dualvals = zeros(nrduals,1);
%     for i=1:nrduals
%         dualvals(i) = value(dual(probability_constraints(i)));
%     end

                
    funcoutput = value(objective);
end

