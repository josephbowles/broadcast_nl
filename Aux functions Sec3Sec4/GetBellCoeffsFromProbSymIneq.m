function [bellcoeffs, constantterm] = GetBellCoeffsFromProbSymIneq(symineq,ins,outs)
% WARNING this only works for 3 parties

constantterm = 0; % b(x,y,z,a,b,c)*p(x,y,z,a,b,c) + constantterm

dims = num2cell([ins,outs]);
bellcoeffs = zeros(dims{:}); % call as bellcoeffs(x,y,z,a,b,c) for 3 parties

nrparties = length(ins);

[C,T] = coeffs(symineq);
nrterms = length(T);

for i=1:nrterms

try 
    double(T(i)); % some terms T(i) are of the form sym('1.1'), they are 
                  % numbers but symbolic. the only way i found to 
                  % recognize them is to convert them to double.
                  % if the variable is sym but not numeric, like p_111_123
                  % then double will give an error and we catch it by doing
                  % the original algorithm, if there is no error then just
                  % assign the coefficient C(i)az to constantterm
    constantterm = constantterm + C(i);
catch
    % assuming input of the form p_abc_xyz and no integer greater than 10
    nameinput = char(T(i));
    nrchars = length(nameinput);

    outputs = zeros(1,nrparties);
    firstidx = 3;
    for idx=1:nrparties
        a = str2double(nameinput(firstidx));
        outputs(idx) = a;
        firstidx = firstidx + 1;
    end
    firstidx = firstidx + 1;
    settings = zeros(1,nrparties);
    for idx=1:nrparties
       x = str2double(nameinput(firstidx));
       settings(idx) = x;
       firstidx = firstidx + 1;
    end

    indexes = num2cell([settings, outputs]);
    bellcoeffs(indexes{:})=C(i);
end % end of try catch

end

end

