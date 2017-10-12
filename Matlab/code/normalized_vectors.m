function d =normalized_vectors(h1, h2)
%   NORMALIZED_VECTORS computes the Normalized Vectors similarity between
%   the two histograms

d= dot(h1./norm(h1),h2./norm(h2));
end