function d = euclideanZ(h1, h2)
%   EuclideanZ compute the euclidean distance after taking the z-score
%   standardizations of h1 and h2
%

h1_=(h1-mean(h1))/std(h1);
h2_=(h2-mean(h2))/std(h2);

d = euclidean(h1_,h2_);
end