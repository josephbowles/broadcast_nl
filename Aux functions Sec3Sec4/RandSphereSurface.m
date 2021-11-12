function vec = RandSphereSurface(dim)
    v = randn(dim,1); 
    n = norm(v);
    if n>1e-6
        vec = v/n;
    else
        vec = RandSphereSurface(dim);
    end
end

