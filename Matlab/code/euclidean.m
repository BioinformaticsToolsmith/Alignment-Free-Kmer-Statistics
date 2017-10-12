function d = euclidean(h1, h2)
%   EUCLIDEAN computes the Euclidean Distance for the two histograms h1 and h2

d = sqrt(sum((h1-h2).^2));
end