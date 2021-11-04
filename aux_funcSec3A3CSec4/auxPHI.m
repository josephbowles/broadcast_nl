function PHI = auxPHI(dim)
    PHI = zeros(dim*dim,1);
    basis=eye(dim);
    for i=1:dim
       PHI = PHI + kron(basis(:,i),basis(:,i));
    end
end