function f=loop_vector(dims) 

%Creates the vector containing all combinations of a_1...a_N, where a_1 is in 1:dims(1), a_2 in 1:dims(2), etc. 

%Useful when one wants to loop over variables a_1...a_N, going from 1:n1, 1:n2, ... ,1:nN

%'dims' is a 1xN vector, containing the n1, n2,...,nN (that is, the number of each 'a_k' one wants to sum over)


N=length(dims);

%Create cell of sets A1...AN: 

C=cell(1,N); 

for k=1:N
    
    C{k}=1:dims(k);
    
end


f=allcomb(C{:});

end