function d = EMD(h1, h2)
%   EMD computes Earth Mover's Distance for the two histograms h1 and h2

d =sum(abs(cumsum(h1)-cumsum(h2)));
end