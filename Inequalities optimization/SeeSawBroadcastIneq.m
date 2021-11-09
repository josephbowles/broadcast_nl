function f = SeeSawBroadcastIneq(Ineq,rho,dims,WhichParties,WhichChois,Choi,Meas,Precision)

%Optimizes broadcast inequality 'Ineq' (for state 'rho' of local dimensions 'dims') over measurements of parties 'WhichParties' and Choi states 'WhichChois' (rest is set from 'Choi' and 'Meas'), with (optional) seesaw precision 'Precision'

%Uses a seesaw of SDPs, the initial point is taken in 'Choi' and 'Meas'


%Inputs:

%'Ineq' is a 1 x oa*ob*oc*...*nx*ny*nz*... vector, interpreted as the coefficients c(abc...|xyz...) = [c(000..|000..) c(00..1|000) ... ]) of the inequality

%'rho' is a (dA*dB*dC..)x(dA*dB*dC..) PSD and trace-1 matrix, representing the (initial) multipartite quantum state

%'dims' is a 1Xn vector, containing the (initial) local dimensions: dims = [dA dB dC ...]

%'WhichParties' is a vector of integer, indicating the (final) parties's labelings one wants to optimize over

%'WhichChoi' is a vector of integer, indicating the (initial) parties's labelings one wants to optimize over

%'Choi' is a 1xn cell, each entry containing a Choi state (that is, a PSD matrix T such that PartialTrace(T,2) = eye(d1) )

%'Meas' is a 1xN cell, each entry containing local POVMs in table form: Meas = { A(dA,dA,oa,nx), B(dB,dB,oa,ny), ... }
%That is, A(dA,dA,oa,nx) is a dA x dA x oa x nx matrix, such that
%A(:,:,a,x)>=0 and sum_n(A(:,:,:,x),3)=eye(dA) for all x=1,..,nx

%Optional input:
%'Precision' is a scalar: the seesaw stops when the difference between two
%succesive values is smaller than 'Precision'

%Outputs:

%'f' is a 1x2 cell, f{1} containing the maximal value of the inputed
%inequality, f{2} the corresponding updated measurements 'Meas' (same
%format as input), f{3} the corresponding updatet Choi states 'Choi' (same
%format as input)



if nargin < 8
    
    Precision = 10^-4;
    
end



%Number of parties

N=length(Meas);


%Repeats the whole seesaw until two succesives value are closer than
%'Precision'

gap = Inf;
CurrentValue=Inf;

while gap > Precision
    
    for k=1:N
        
        
        if ismember(k,WhichChois)
            
            g = MaxIneqBroadcastChoi(Ineq,rho,dims,k,Choi,Meas);
            
            Choi{k}=g{2};
            
        end
        
        if ismember(k,WhichParties)
            
            g = MaxIneqBroadcastMeas(Ineq,rho,dims,k,Choi,Meas);
            
            Meas{k}=g{2};
            
        end
        
        
        
    end
    
    
    
    gap = abs(CurrentValue - g{1}) ;
    
    
    CurrentValue = g{1};
    
end




f = {CurrentValue,Meas,Choi};



end