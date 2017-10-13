function d = hellinger( h1, h2 )
%   HELLINGER computes the hellinger distance between the two histograms

h1_ave=mean(h1);
h2_ave=mean(h2);

d = sqrt(2*sum((sqrt(h1./h1_ave)-sqrt(h2./h2_ave)).^2));
end