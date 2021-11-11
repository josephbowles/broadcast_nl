function subscript_vector = ind2subv(siz, linear_index)
%IND2SUBV return vector position of a linear index into an ND-array
%
% subscript_vector = ind2subv(siz, linear_index)
%
% Inputs:
%                   siz Dx1 or 1xD   size of array to index
%          linear_index Nx1 or 1xN   required elements using Fortran order
%
% Outputs:
%     subscript_vector      NxD      required elements as subscript row-vectors
%
% Matlab's ind2sub turns linear scalar indexes into the subscripts corresponding
% to the row, col, 3rd-dim, ... of the element. The subscript for each dim
% is returned as a separate argument. This version returns a single argument:
% each row is a vector containing all the subscripts for the corresponding
% linear_index.
%
% See also: subv2ind, ind2sub, sub2ind

% Iain Murray, November 2007, March 2011
% This technique appeared earlier in the following forum post:
%
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=5476&objectType=File
% Date: 2004-10-04
% From: Mukhtar Ullah (mukhtar.ullah@informatik.uni-rostock.de)

[out{1:length(siz)}] = ind2sub(siz, linear_index(:));
subscript_vector = cell2mat(out);
end