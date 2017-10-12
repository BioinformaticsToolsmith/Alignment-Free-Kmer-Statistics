function d = jensen_shannon( h1 , h2 )
%   JENSEN_SHANNON computes the Jensen-Shannon divergence between the two
%   histograms

p1=h1./sum(h1);
p2=h2./sum(h2);
m=p1+p2;
%Works better with probabilities instead of frequencies of k-mers
d =(1/2)*sum(p1.*log(2*p1./m)+p2.*log(2*p2./m));

%histogram with itself yields infinity?
end