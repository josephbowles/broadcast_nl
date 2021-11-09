function rho_out=ApplyChoiMap(rho,T,DIM,position)
%T is the Choi operator of the target map
%DIM is a vector with the dimensions of all systems
%position is the position which you desire to perform T.

%Let illustrate the problem with tripartite states with spaces labelled by
%1, 2 and 3, and position=3
%The output of this algorithm will be
%\tr_{1'2'3'} (|Id><Id|_{11'} \otimes |Id><Id|_{22'} \otimes T_{33'}
%  transpose (rho_{123})  \otimes Id_{1'2'3'}
%Where the 1',2',3' are the spaces are auxiliary spaces
%Since matlab does not accept 1',2', and 3',
% we say that 1' is 4, 2' is 5, and 3' is 6


k=size(DIM,2); %find the number of parties
d_out=size(T,1)/DIM(position); %Evaluate the output dimension of T
if k==2
    switch position
        case 1
            rho=PermuteSystems(rho,[2 1],DIM);
rho_out = PartialTrace(kron(PartialTranspose(rho,2,[DIM(2) DIM(1)]),eye(d_out))*kron(eye(DIM(2)),T),2,[DIM(2) DIM(1) d_out]);
            rho_out=PermuteSystems(rho_out,[2 1],[DIM(2) d_out]);
        case 2
rho_out = PartialTrace(kron(PartialTranspose(rho,2,DIM),eye(d_out))*kron(eye(DIM(1)),T),2,[DIM d_out]);
	
    end
    
else




%Construct the total Choi operator which will be performed in rho
%Choi total will respect the order 11' 22' 33' 44'
%Construct the DIM_total in order 11' 22' 33' 44'
%Construct rho_total in order 1234  1'2'3'4'
%Construct the DIM_total_rho in order  1234  1'2'3'4'
ChoiTotal=1;
rho_total=rho;
DIM_total=0;
DIM_total_rho=DIM;
for i=1:k
    if i==position
        ChoiTotal=kron(ChoiTotal,T);
        rho_total=kron(rho_total,eye(d_out));
        DIM_total=[DIM_total DIM(i) d_out];
        DIM_total_rho=[DIM_total_rho d_out];
    else
        ChoiTotal=kron(ChoiTotal,IsotropicState(DIM(i),1))*DIM(i);
        rho_total=kron(rho_total,eye(DIM(i)));
        DIM_total=[DIM_total DIM(i) DIM(i)];
        DIM_total_rho=[DIM_total_rho DIM(i)];
    end
end
DIM_total=DIM_total(2:end); %Remove the initial 0 from DIM_total


permute_vector4rho=0; %Construct the permutation vector to write rho as 11' 22' 33' 44' 
%instead of 11' 22' 33' 44'
for i=1:k
    permute_vector4rho=[permute_vector4rho i k+i]; 
end
permute_vector4rho=permute_vector4rho(2:end); %Remove the initial 0 permute_vector

rho_total_permuted = PermuteSystems(rho_total,permute_vector4rho,DIM_total_rho); %Construct rho in the order 11'22'33'44'

aux_systems=2*(1:k)-1; %Creates a vector that stocks all auxiliary systems

rho_out=PartialTrace(ChoiTotal*transpose(rho_total_permuted),aux_systems,DIM_total);

end
end
