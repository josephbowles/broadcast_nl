function f=Choijoe(U,d1) 

%Creates the Choi state corresponding to the map that sends the state of Bob0 to Bob1 (both of dimension 'd1'), while Bob2 has a ket0, and applies a unitary 'U' between Bob1 and Bob2

%'U' is a d1*d2 x d1*d2 unitary matrix (that is, one needs size(U,1)/d1 to be an integer)
%'d1' is an integer 


d2=size(U,1)/d1;

K=zeros(d1*d2,d1); 

for k=1:d1
    
    K((k-1)*d2+1,k)=1;
    
end



f=ChoiMatrix({U*K}); 

% %check
% TOT=kron(transpose(rho),eye(d1*d2))*T ; 
% 
% rhot=PartialTrace(TOT,1,[d1 d1*d2]);
% 
% norm(rhot-U*K*rho*K'*U')


end
