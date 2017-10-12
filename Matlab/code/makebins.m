function [T] = makebins(Set,Y,bins)
%   MAKEBINS divides a given data set into cell bins
%   -Set is the Data set (usually consisting of alignment identity score
%   values and their respective sequences
%   -Y is the result of running MATLAB's discretize function on Set
%   -bins is the number of bins that the data will be divided into

T=cell(bins,1);
for i=1:size(Y,1)
   num=Y(i,1);
   elt=Set(i,:);
   T{num,1}=[T{num,1}; elt];
end


end
