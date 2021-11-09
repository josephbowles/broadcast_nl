function f=LHSbroadcast_3BobsCG(sigrho,sigN,scenario,DB_CG)

%Finds the largest q such that q*'sigrho'+(1-q)*'sigN' is unsteerable in the broadcast scenario (with 3 Bobs), characterized by 'scenario', and using vertices 'DB' for the  Bobs

%Inputs:

%'sigrho' is a dA x dA x ob x oc x od x ny x nz x nw matrix, interpreted as a broadcast assemblage
%'sigN' is a dA x dA x ob x oc x od x ny x nz x nw matrix, interpreted as a broadcast assemblage
%'scenario' is of the form [ob oc od |ny nz nw]
%'DB' is a L x N_CG matrix, each line of which is interpreted as a distribution DB(k,:) = p(b1b2|y1y2) = [p(00|00) p(01|00) ... p(00|01) ..], using Collins-Gisin notation (N_CG is the Collins-Gisin dimension)


%Outputs:

%'O' is the value is the solution of the SDP (that is, the value of the objective, the best visibility 'q')
%'F' is the dual variables, that is, the optimal steering inequality for
%assemblages 'sigrho', 'sigN' (will be given in Collins-Gisin notation)



%Scenario

ob=scenario(1); ny=scenario(4);

oc=scenario(2); nz=scenario(5);

od=scenario(3); nw=scenario(6);


dA=size(sigrho,1);



%convert vertices to a table
DB=DB_CG;

L=size(DB,1);

Df=zeros(ny,nz,nw,ob,oc,od,L);

for k=1:L
    
    Df(:,:,:,:,k)=convP2Table(DB(k,:),[ny nz nw ob oc od]);
    
end



%%%%SDP

cvx_begin quiet

%variables:

variable q nonnegative
variable X(dA,dA,L) hermitian semidefinite;
dual variable F

% for k=1:L
%
%     sig(:,:,k) == semidefinite(dA,dA);
%
% end




q <= 1;


%LHS decomposition

%Full 
S=cvx(zeros(dA,dA,ob,oc,od,ny,nz,nw));
for y=1:ny
    for z=1:nz
        for w=1:nw
            
            for b=1:ob-1
                for c=1:oc-1
                    for d=1:od-1
                        
                        
                        for k=1:L
                            S(:,:,b,c,d,y,z,w) = S(:,:,b,c,d,y,z,w) + X(:,:,k)*Df(b,c,d,y,z,w,k);
                        end
                        
                    end
                end
            end
            
            
        end
    end
end



%Marginals


%
F: S == q*sigrho + (1-q)*sigN  ;




%objective function

O = q;



maximize(O);

cvx_end



f={O,F};


end