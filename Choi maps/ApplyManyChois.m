function rho_out=ApplyManyChois(rho,Choi,dims)

%Applies all Choi maps in 'Choi' to state 'rho', of local dimension 'dims'

%Inputs:
%'rho' is a (dA*dB*dC..)x(dA*dB*dC..) PSD and trace-1 matrix, representing a multipartite quantum state

%'dims' is a 1Xn vector, containing the local dimensions: dims = [dA dB dC ...]

%'Choi' is a 1xn cell, each entry containing a Choi state (that is, a (d1*d2)x(d1*d2) PSD matrix T such that PartialTrace(T,[d1 d2],2) = eye(d1) )


%Outputs:
%'rho_out' is the final state



%final state
dimsc=dims;
for k=1:length(dims)
    rho=ApplyChoiMap(rho,Choi{k},dimsc,k);
    dimsc(k)=dims(k)*size(rho,1)/prod(dimsc); 
end


rho_out = rho;

end