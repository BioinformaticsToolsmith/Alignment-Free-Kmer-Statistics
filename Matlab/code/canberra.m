function d =canberra(h1, h2)
%   CANBERRA computes the canberra distance for two histograms h1 and h2

d=sum(abs(h1-h2)./(h1+h2));
end