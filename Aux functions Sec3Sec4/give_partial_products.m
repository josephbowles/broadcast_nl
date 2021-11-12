function [partial_products] = give_partial_products(povms, bellcoeffs, party, ins, outs)

nrparties = length(ins);

%for party_to_ignore=1:nrparties
party_to_ignore = party;
    partial_products = cell(ins(party_to_ignore), outs(party_to_ignore));
    allbutone = logical(ones(nrparties,1)); 
    allbutone(party_to_ignore) = false;
    ins_without_party_to_ignore = ins(allbutone);
    outs_without_party_to_ignore = outs(allbutone);
    
    dimen = prod(outs_without_party_to_ignore);
    for x=1:ins(party_to_ignore)
        for a=1:outs(party_to_ignore)
            partial_products{x,a} = 0;
            for y=1:ins_without_party_to_ignore(1)
                for z=1:ins_without_party_to_ignore(2)
                    for b=1:outs_without_party_to_ignore(1)
                        for c=1:outs_without_party_to_ignore(2)
                            if party_to_ignore == 1
                                term = bellcoeffs(x,y,z,a,b,c) * kron(povms{2}{y}{b},povms{3}{z}{c});
                            elseif party_to_ignore == 2
                                term = bellcoeffs(y,x,z,b,a,c) * kron(povms{1}{y}{b},povms{3}{z}{c});
                            elseif party_to_ignore == 3
                                term = bellcoeffs(y,z,x,b,c,a) * kron(povms{1}{y}{b},povms{2}{z}{c});
                            else
                               error('more than 3 parties not supported'); 
                            end
                            partial_products{x,a} = partial_products{x,a} + term;
                        end
                    end
                end
            end
        end
    end
    
%end

end