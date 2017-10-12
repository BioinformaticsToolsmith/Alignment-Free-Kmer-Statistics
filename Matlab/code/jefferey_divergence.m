function d = jefferey_divergence ( h1 , h2 )
%   JEFFEREY_DIVERGENCE computes the jefferey divergence between the two
%   histograms

p1=h1./sum(h1);
p2=h2./sum(h2);
%Works better with probabilities instead of frequencies of k-mers
d =sum((p1-p2).*log(p1./p2));

end