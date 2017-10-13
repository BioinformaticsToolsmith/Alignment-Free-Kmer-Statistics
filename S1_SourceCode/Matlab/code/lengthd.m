function d = lengthd(h1,h2)
% LENGTHD computes the length difference between the two histograms
%
%

l1=sum(h1);
l2=sum(h2);

d=abs(l1-l2);

end