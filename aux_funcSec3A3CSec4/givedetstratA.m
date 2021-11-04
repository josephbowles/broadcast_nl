function det = givedetstratA(out,in)
    nrOfOuts = out;
    nrOfIns = in;
    detpoints = nrOfOuts^nrOfIns;
    
    % combinations with arbitrary nr of length
    % https://es.mathworks.com/matlabcentral/answers/433424-how-to-get-the-combinations-of-elements-of-two-arrays
    C = cell(1,nrOfIns);
    for i = 1:nrOfIns
       C{i} = 1:out; 
    end
    D = C;
    [D{:}] = ndgrid(C{:});
    Z = cell2mat(cellfun(@(m)m(:),D,'uni',0));
    
    % AlternaCtive way to do the above:
    %outs=[2,2,2,2];
    %cartproductOUT = ind2subv(outs, 1:prod(outs(:)));
    
    det = zeros(detpoints,nrOfIns,nrOfOuts);
    
    lam = 1;
    for i1=1:(length(Z(:,1)))
        atuple = Z(i1,:);
        for x = 1:nrOfIns
            det(lam,x,atuple(x))=1;
        end
        lam = lam + 1;
    end
end
