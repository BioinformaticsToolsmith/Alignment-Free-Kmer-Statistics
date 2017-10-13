function d = markov(h1, h2)
%   MARKOV calculates the MARKOV similarity between the two histograms
%   note: Markov is not symmetric, one sequence becomes the model
%   -h1 is the sequence
%   -h2 is the model
    cond=log(build_cond(h2));
    d=sum((h1-1).*cond);
end

