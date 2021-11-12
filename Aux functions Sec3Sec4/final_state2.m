function state = final_state2(inistate,channel_left, channel_right)
    % This is a modification of the unnumbered function for the 4 party scenario
    % with 2 alices and 2 bobs.

    % TODO put dimA, dimB, dimB1, dimB2 as function inputss
    
    dimA=2;
    dimB=2;
    dimB1=2;
    dimB2=2;
    
       
    reshaped_state = reshape(inistate,dimA,dimB,dimA,dimB);
    state = 0;
    for i=1:dimA
        for j=1:dimB
            for k=1:dimA
                for l=1:dimB
                    state = state + reshaped_state(i,j,k,l) * ...
                                    kron(  ApplyMap(ketbra(i,k,dimA),ChoiMatrix(channel_left)), ...
                                           ApplyMap(ketbra(j,l,dimB),ChoiMatrix(channel_right)) );
                end 
            end
        end
    end
       
    
end


