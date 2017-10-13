function d= SimMM(h1, h2)
%   SIMMM computes the Similarity of Markov models for the two histograms

d=1-exp((d_markov(h1,h2)+d_markov(h2,h1))/2);

end


function d = d_markov(h1,h2)
%   D_MARKOV is a helper function that makes use of markov

d=(1/sum(h2-1))*log((markov(h2,h1))/markov(h2,h2));

end

