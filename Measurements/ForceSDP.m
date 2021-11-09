function f=ForceSDP(A)

%"Makes" 'A' SDP

%First define A2=(A+A')/2, then remove the subspace corresponding to negative eigenvalues of A2


A=(A+A')/2; 

[V,D] = eig(A); 

D(D < 0) = 0; 


f=V*D*V';

end