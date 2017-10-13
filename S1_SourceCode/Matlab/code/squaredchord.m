function d =squaredchord(h1, h2)
%   SQUARED_CHORD computes the squared chord distance between the two
%   histograms

d=sum((sqrt(h1)-sqrt(h2)).^2);
end