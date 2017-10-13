function d = harmonic_mean( h1 , h2 )
% HARMONIC_MEAN computes the harmonic mean similiarity between the two
% histograms
%
%

p1=h1./sum(h1);
p2=h2./sum(h2);
%Works better with probabilities instead of frequencies of k-mers
d  = 2.*sum((p1.*p2)./(p1+p2));
end