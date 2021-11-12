function [belloperator] = give_Bell_operator(bellcoeffs, povms, ins, outs)
belloperator = zeros(prod(outs(:)),prod(outs(:)));
for x=1:ins(1)
    for y=1:ins(2)
        for z=1:ins(3)
            for a=1:outs(1)
                for b=1:outs(2)
                    for c=1:outs(3)
                        term = Tensor(povms{1}{x}{a}, ...
                                        povms{2}{y}{b},...
                                        povms{3}{z}{c});
                        belloperator = belloperator + bellcoeffs(x,y,z,a,b,c)*term;
                    end
                end
            end
        end
    end        
end
end