function [ia_povms] = give_ia_povms(POVMS, ins, outs)

nrparties = length(ins);

ia_povms = cell(nrparties, max(ins), max(outs)); 
for party=1:nrparties
    for in=1:ins(party)
        for out=1:outs(party)
            ia_povms{party,in,out} = POVMS{party}{in}{out}(find(triu(ones(outs(party)))));
        end
    end
end
end