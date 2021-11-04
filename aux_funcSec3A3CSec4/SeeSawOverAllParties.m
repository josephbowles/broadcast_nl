function [finalPproj,finalobj,finalchannel] = SeeSawOverAllParties(bellcoeffs, initialState, initialPovms, initialChannel)

    ABS_TOL   = 1e-6;
    MAX_ITER  = 200;
    abschange = 1e6;
    objval    = -1e6;

    povms   = initialPovms;
    channel = initialChannel;

    assert(mod(length(size(bellcoeffs)),2)==0,"There should be as many inputs as outputs.");
    nrparties = length(size(bellcoeffs))/2;
    
    partyidx = floor(1+4*rand(1)); % start randomly
    iteration = 1;
    while abschange>ABS_TOL && iteration <= MAX_ITER
        if partyidx < nrparties+1
            output_state = final_state(initialState, channel);
            [newPproj,newobjval,problemStatus] = SeeSawOverASingleParty(partyidx, output_state, bellcoeffs, povms);
            povms = newPproj;
            
            abschange = abs(newobjval - objval);
            
            %fprintf("iter=%d, partyidx=%d, bellineqvalue=%f, abschangeobjval=%g\n", iteration, partyidx, newobjval, abschange);
            objval = newobjval;
        else
            % now we do the channel
            [newchannel,newobjval,problemStatus] = SeeSawOverChannel(initialState, bellcoeffs, povms);
            channel = newchannel;
                        
            abschange = abs(newobjval-objval);
            %fprintf("iter=%d, channelq=%c, bellineqvalue=%f, abschangeobjval=%g\n", iteration, "~", newobjval, abschange);
   
            objval = newobjval;
        end
        partyidx = mod(partyidx,nrparties+1)+1;
        iteration = iteration + 1;
    end
    finalPproj = povms;
    finalobj = objval;
    finalchannel = channel;
end