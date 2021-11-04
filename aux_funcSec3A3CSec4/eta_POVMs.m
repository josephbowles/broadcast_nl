function newPOVMs= eta_POVMs(povms, party, strat, eta, ins, outs)
% Auxiliary function used for calculating the noisy distribution resuling
% from assuming detectors sometimes fail.


%assert(0 <= eta <= 1, "Invalid value for eta.");

newPOVMs = povms;

ze = zeros(outs(party));
id = eye(outs(party));

if party==1
    if strat==1
        det = {{id,ze},{id,ze},{id,ze}};
    elseif strat==2
        det = {{id,ze},{id,ze},{ze,id}};
    elseif strat==3
        det = {{id,ze},{ze,id},{id,ze}};
    elseif strat==4
        det = {{id,ze},{ze,id},{ze,id}};
    elseif strat==5
        det = {{ze,id},{id,ze},{id,ze}};
    elseif strat==6
        det = {{ze,id},{id,ze},{ze,id}};
    elseif strat==7
        det = {{ze,id},{ze,id},{id,ze}};
    elseif strat==8
        det = {{ze,id},{ze,id},{ze,id}};
    else
        error("invalid strat");
    end
else
    if strat==1
        det = {{id,ze},{id,ze}}; % x=0->a=0, x=1->a=0
    elseif strat==2
        det = {{id,ze},{ze,id}}; % x=0->a=0, x=1->a=1
    elseif strat==3
        det = {{ze,id},{id,ze}}; % x=0->a=1, x=1->a=0
    elseif strat==4
        det = {{ze,id},{ze,id}}; % x=0->a=1, x=1->a=1
    else
        error("invalid strategy");
    end
end
    
for x = 1:ins(party)
   for a = 1:outs(party)
       newPOVMs{party}{x}{a} = eta * povms{party}{x}{a} + (1-eta)* det{x}{a};
   end
end

end

