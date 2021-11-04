function out = BlochComponents(obs)
    sig0 = eye(2);
    sig1 = [0+0i,1+0i;1+0i,0+0i];
    sig2 = [0+0i,-1i;1i,0+0i];
    sig3 = [1+0i,0+0i;0+0i,-1+0i] ;
    
    r0 = trace(obs' * sig0);
    r1 = trace(obs' * sig1);
    r2 = trace(obs' * sig2);
    r3 = trace(obs' * sig3);
    
    out = [r0, r1, r2, r3]/2;
end

