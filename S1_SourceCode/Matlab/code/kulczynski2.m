function d =kulczynski2(h1, h2)
%   KULCZYNSKI2 computes Kulczynski Similarity 2 between the two
%   histograms

h1_ave=mean(h1);
h2_ave=mean(h2);

d=(size(h1,2)/2)*((1/h1_ave)+(1/h2_ave))*sum(min(h1,h2));
end