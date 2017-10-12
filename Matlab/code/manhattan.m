function d = manhattan(h1, h2)
% MANHATTAN Manhattan distance between histograms
d =sum(abs(h1-h2));
end