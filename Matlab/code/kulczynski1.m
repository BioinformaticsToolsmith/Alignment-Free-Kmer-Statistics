function d =kulczynski1(h1, h2)
%   KULCZYNSKI1 computes Kulczynski Similarity 1 between the two
%   histograms

d=sum(abs(h1-h2)./(min(h1,h2)));
end