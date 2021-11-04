function out = checkPOVMsAreGood(povms, ins, outs, tol)

for p = 1:length(ins)
   for x=1:ins(p)
       summ = 0;
      for a=1:outs(p)
          summ = summ + povms{p}{x}{a};
          if ~IsPSD(povms{p}{x}{a},tol)
              warning("POVM not positive! min(eig)=%g\n", min(eig(povms{p}{x}{a})));
          end
      end
      if norm(summ-eye(2),'fro')>tol
          warning("POVMs don't summ up to 1! dist_to_eye(2)=%g\n", norm(summ-eye(2)));
      end
   end
end

out = true;
end

