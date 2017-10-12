function d = rre_k_r(h1, h2)
%   RRE_K_R computes the Revised Relative Entropy similiarty for the two
%   histograms

prob1=build_cond(h1);
prob2=build_cond(h2);
prob_sum=prob1+prob2;

d1= sum(prob1.*log(2*prob1./(prob_sum)));
d2=sum(prob2.*log(2*prob2./(prob_sum)));
d=(d1+d2)/2;



end