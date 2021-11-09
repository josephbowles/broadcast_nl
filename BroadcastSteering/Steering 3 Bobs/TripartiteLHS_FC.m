function f=TripartiteLHS_FC(sigrho,sigN,scenario,Vertices)

%Using the full-correlation notation, finds the largest q such that q*'sigrho'+(1-q)*'sigN' is unsteerable in the broadcast scenario, characterized by 'scenario', and using vertices 'Vertices' for Bob, Charlie and David


%'sigrho' is a dA x dA x (ny+1) x (nz+1) x (nw+1) matrix, interpreted as an assemblage in the broadcast scenario

%'sigN' is a dA x dA x (ny+1) x (nz+1) x (nw+1) matrix, interpreted as an assemblage in the broadcast scenario

%'scenario' is of the form [ny nz nw]

%'Vertices' is a (ny+1) x (nz+1) x (nw+1) x L matrix, corresponding to
%vertices given in the full-correlator form, that is, Vertices(y,z,w,k) is
%the output for vertex k and choice of inputs y,z and w, and y+1 means the
%idendity measurement on y, etc.


%Extracts scenario 

ny=scenario(1); nz=scenario(2); nw=scenario(3); 

dA=size(sigrho,1);


L = size(Vertices,4);


%%%%SDP

%cvx_begin quiet

cvx_begin



%variables:

variable q nonnegative
variable X(dA,dA,L) hermitian semidefinite;
dual variable F


q <= 1;


%LHS decomposition

S=cvx(zeros(dA,dA,ny+1,nz+1,nw+1));

for y=1:ny+1
    for z=1:nz+1
        for w=1:nw+1
            
            
            
            for k=1:L
                
                S(:,:,y,z,w) = S(:,:,y,z,w) + X(:,:,k)*Vertices(y,z,w,k);
                
            end
            
            
        end
    end
end




F: S == q*sigrho + (1-q)*sigN  ;



%objective function

O = q;



maximize(O);

cvx_end



f={O,F};


end