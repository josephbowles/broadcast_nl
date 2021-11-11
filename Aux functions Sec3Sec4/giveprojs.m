function y = giveprojs(op)
    [vec,val] = eig(op);
    s = size(val);
    nrvecs = s(1);
    y = cell(nrvecs);
    for i=1:nrvecs
        v = vec(:,i);
        y{i} = v*v'/(norm(v)^2);
    end
end