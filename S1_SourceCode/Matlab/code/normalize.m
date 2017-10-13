function feature = normalize(f)
%   NORMALIZE scales each of the statistics between 0 and 1

smallest=min(f);

feature= (f-smallest)/(max(f)-smallest);

end

