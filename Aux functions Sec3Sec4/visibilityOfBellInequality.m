function [visibility, LPstatus] = visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform)

%% LP way
% assert(mod(length(size(bellcoeffs)),2)==0,"There should be as many inputs as outputs.");
% dims = size(bellcoeffs);
% nrparties = length(dims)/2;
% ins = dims(1:nrparties);
% outs = dims(nrparties+1:end);
% 
% alpha = sdpvar(1);
% vis_constraints = [alpha>=0, alpha<=1];
% 
% p_noisy = (1-alpha)*p_entangled + alpha*p_uniform;
% % objective = sum(bellcoeffs .* p_noisy, size(bellcoeffs));
% % aux = bellcoeffs .* p_noisy;
% % auxsize=size(aux);
% % while prod(auxsize(:)) ~= 1
% %    aux = sum(aux); 
% % end
% 
% objective = 0;
% coords = ind2subv(size(bellcoeffs), 1:prod(dims(:)));
% for idx = 1:size(coords,1)
%     coords_choice = num2cell(coords(idx,:));
%     objective = objective + bellcoeffs(coords_choice{:}) * p_noisy(coords_choice{:});
% end
% %sdisplay(objective);
% constraints = [vis_constraints, objective <= localbound];
% 
% optsol = optimize(constraints, -objective, sdpsettings('solver','mosek','verbose', 0));
% LPstatus = optsol.problem;
% if optsol.problem ~= 0
%     warning('Infeasibility in inequality visibility LP.');
%    %error('check why this LP is not working'); 
% end
% 
% visibility = value(alpha);

%% Manual way

aux1 = bellcoeffs .* p_entangled;
c1 = sum(aux1(:));

aux2 = bellcoeffs .* p_uniform;
c2 = sum(aux2(:));

comparison_tolerance = 1e-8;
if abs(c1-c2) < comparison_tolerance
    visibility = 0;
    LPstatus = 1;
    diff_p1p2 = norm(p_entangled(:) - p_uniform(:));
    if diff_p1p2 < comparison_tolerance
        fprintf("! The probability distributions are too close. tol: %f\n", comparison_tolerance);
    else
        fprintf("! The Bell values of the probability distributions are too close. tol: %f\n", comparison_tolerance);
    end
    return;
end

if c1 <= localbound && c2 <= localbound
    visibility = 0;
    % technically it is an infeasible problem, but for this case it's best 
    % not to raise this flag as this is a common situation. 
    % it usually gives visibility=1 or somthing like that
    LPstatus = 0; 
    
    %warning("Neither p1 nor p2 violate the Bell inequality. Returning 0.");
    return;
elseif c1 >= localbound && c2 >= localbound
    visibility = 0;
    LPstatus = 1;
    warning("Both p1 and p2 violate the Bell inequality. Returning 0.\n");
    return;
elseif c1 <= localbound && c2 >= localbound
    visibility = 0;
    LPstatus = 1;
    warning("The code convention is that p1 is outside the local set and p2 inside. You probably didn't optimize over channels and measurements. Returning 0.\n");
    return;
else
   % only thing left is c1 >= localbound and c2 <= localbound
   visibility = (c1-localbound)/(c1-c2);
   LPstatus = 0;
end

end



