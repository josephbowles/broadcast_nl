function state = NoisyPartiallyEntangled(p,xi,identity_placement)
    %e1 = [1;0];
    %e2 = [0;1];

    %e1e1 = Tensor(e1,e1);
    %e1e2 = Tensor(e1,e2);
    %e2e1 = Tensor(e2,e1);
    %e2e2 = Tensor(e2,e2);
    
    e1e1 = [1;0;0;0];
    e2e2 = [0;0;0;1];
    
    psi_xi = cos(xi) * e1e1 + sin(xi) * e2e2;
    rho_xi = psi_xi * psi_xi'/ (psi_xi' * psi_xi);
    
    if identity_placement == 'A'
        rhoB = PartialTrace(rho_xi, 1, [2,2]);
        state = (1-p) * rho_xi + p * kron(eye(2)/2, rhoB);
    elseif identity_placement == 'B'
        rhoA = PartialTrace(rho_xi, 2, [2,2]);
        state = (1-p) * rho_xi + p * kron(rhoA, eye(2)/2);
    else
       error("Invalid 'identity_placement' input."); 
    end
end

