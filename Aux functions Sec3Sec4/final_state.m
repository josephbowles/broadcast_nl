function state = final_state(inistate,channel)
    % TODO put dimA, dimB, dimB1, dimB2 as function inputss
    dimA=2;
    dimB=2;
    dimB1=2;
    dimB2=2;
    
    %bigchoi = ChoiMatrix({kron(eye(2),channel)});
    
    biggerstate = kron( inistate.', eye(dimA*dimB1*dimB2) );
    
    Phi = auxPHI(dimA);
    biggerchannel = kron(Phi*Phi',ChoiMatrix(channel));
    swapop=Tensor(eye(dimA),SwapOperator(2),eye(dimB1*dimB2));
    biggerchannel = swapop * biggerchannel * swapop';
    
    %using auxSwap DOESNT WORK
    %biggerchannel = auxSwap(biggerchannel, [dimA,dimA,dimB,dimB1,dimB2], [2,3]);
    
    state = PartialTrace(biggerchannel*biggerstate, [1,2], [dimA,dimB,dimA,dimB1,dimB2]);
%     PartialTrace(biggerchannel, [3,4,5], [dimA,dimB,dimA,dimB1,dimB2])
%     
%     PartialTrace(bigchoi*biggerstate, [1,2], [dimA,dimB,dimA,dimB1,dimB2])
%     PartialTrace(bigchoi, [3,4,5], [dimA,dimB,dimA,dimB1,dimB2])
%     
%     bool1 = bigchoi == biggerchannel;
%     fprintf("is chanel ok? %d", prod(bool1(:)));
%     
%     out1 = PartialTrace(biggerstate*biggerchannel, [1,2], [dimA,dimB,dimA,dimB1,dimB2]);
%     out2 = ApplyMap(inistate, biggerchannel);
% 
%     
%     reshaped_state = reshape(inistate,dimA,dimB,dimA,dimB);
%     state = 0;
%     for i=1:dimA
%         for j=1:dimB
%             for k=1:dimA
%                 for l=1:dimB
%                     state = state + reshaped_state(i,j,k,l) * ...
%                                     kron( ketbra(i,k,dimA), ApplyMap(ketbra(j,l,dimB),choi));
%                 end 
%             end
%         end
%     end
%     
%     state2=ApplyMap(inistate, bigchoi);
    
end


