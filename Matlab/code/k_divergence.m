function d = k_divergence ( h1 , h2 )
%   K_DIVERGENCE computes the K divergence between the two
%   histograms

p1=h1./sum(h1);
p2=h2./sum(h2);

%Works better with probabilities instead of frequencies of k-mers
d =sum(p1.*log((2.*p1)./(p1+p2)));


end