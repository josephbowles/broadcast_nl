function projs = givePprojRANDgeneral(ins)
    sig1 = [0+0i,1+0i;1+0i,0+0i];
    sig2 = [0+0i,-1i;1i,0+0i];
    sig3 = [1+0i,0+0i;0+0i,-1+0i] ;

    projs={{0}};
    for party = 1:length(ins)
       for pIn = 1:ins(party)
           nP = RandSphereSurface(3);
           P = nP(1) * sig1 + nP(2) * sig2 + nP(3) * sig3;
           projs{party}{pIn} = giveprojs(P); % built-in with outs=[2,2,2], only qubits
       end
    end
end