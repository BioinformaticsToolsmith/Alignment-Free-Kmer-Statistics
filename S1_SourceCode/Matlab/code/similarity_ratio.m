function d =similarity_ratio(h1, h2)
%   SIMILARITY_RATIO computes the Similarity Ratio for the given histograms

product=dot(h1,h2);

d=product/(product+norm(h1-h2));

end