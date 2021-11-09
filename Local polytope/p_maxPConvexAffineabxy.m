function [f,x,exitflag,output,W]=p_maxPConvexAffineabxy(pT,pN,Dc,DA,N_Outputs)

%Finds the largest v such that v*'pT'+(1-v)*'pN' admits a mixed convex-affine
%model in scenario characterized by 'N_Outputs', using convex polytope given by vertices 'Dc', and affine space
%characterized by vertices 'DA'

%Here we assume normalized 'pT','pN' and vertices in 'Dc', 'DA' 

%If 'DA' is the set of local vertices, then basically finds the best model
%over all no-signaling strategies for corresponding group of parties  (see e.g.
%Theorem 1 of https://arxiv.org/pdf/0804.4859.pdf)

%Inputs:
%'pT' is a d x 1 vector, each line of which is interpreted as a behaviour p(ab..|xy..). That is, pT = [p(00..|00..) p(00..1|00..)   ...   p(00...|00...1) p(00..1|0..1) ... ]
%'pN' is a d x 1 vector, each line of which is interpreted as a behaviour p(ab..|xy..). That is, pT = [p(00..|00..) p(00..1|00..)   ...   p(00...|00...1) p(00..1|0..1) ... ]
%'Dc' is a N x d1 matrix, interpreted as a set of N vertices of dimension d1
%'DA' is a L x d2 matrix, interpreted as a set of L vertices of dimension d2
%'N_Outputs' is a 1 x 2 vector, containing the (global) number of outputs of
%both parties. That is, N_Outputs(1) is the number of overall outputs of parties whose vertices are 'Dc', N_Outputs(2) is the number of overall outputs of parties whose vertices are 'DA'

%One must have d1*d2 = d

%Outputs:
%'f' is a scalar, corresponding to the largest number s.t. f*'pT'+(1-f)*'pN' admits a mixed convex-affine
%model



%number of vertices
N=size(Dc,1);
L=size(DA,1);

%local dimension of Party2 
dAff = size(DA,2); 



%overall vertices
V=Concatenate_Vertices(Dc,DA,N_Outputs(1),N_Outputs(2));

%Number of variables
D=L*N+1;


%objective function
F=zeros(D,1);

F(1)=-1;


%equality constraints (v*'pT'+(1-v)*'pN' has to be written as a linear combination of the vertices 'V' )
B=[pN-pT,V'];

c=pN;


%constraint on the positivity of Party2's local strategy 

A = - [zeros(dAff*N,1) kron(eye(N),DA')]; 

%linprog options
% options = optimoptions('linprog','Algorithm','interior-point','Display','off');
% [x,fval,exitflag,output,W]=linprog(F,A,zeros(dAff*N,1),B,c,[],[],[],options);

[x,fval,exitflag,output,W]=linprog(F,A,zeros(dAff*N,1),B,c);


%if there is a problem displays the flag
if (exitflag ~= 1)
    
    f=[-fval,exitflag];
    
else

    f=-fval;

end



end