function d = mismatch( h1, h2 )
%   MISMATCH computes the mismatch distance between the two histograms

d = sum ( h1 ~= h2 );
end