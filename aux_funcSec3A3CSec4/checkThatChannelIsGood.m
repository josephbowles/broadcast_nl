function result = checkThatChannelIsGood(channel, dimx, dimy, tol)

choi = ChoiMatrix(channel);

tracedoutput = PartialTrace(choi, [2], [dimx, dimy]);
diff = tracedoutput - eye(dimx);
if ~(norm(diff) <= tol)
    warning('Choi matrix condition not satisfied. norm(diff)=%g, tol=%g', norm(diff), tol)
end
choieigs = eig(choi);
if ~(all(choieigs >= -tol))
    warning('Choi matrix not semi-definite positive. min(eig(choi))=%g tol=%g', min(choieigs), tol);
end
result=true;
end