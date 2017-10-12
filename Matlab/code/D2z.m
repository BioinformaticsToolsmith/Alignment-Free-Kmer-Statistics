function d = D2z(h1, h2)
% D2Z computes the D2z distance for the two histograms h1 and h2

h1_=(h1-mean(h1))/std(h1);
h2_=(h2-mean(h2))/std(h2);
d=dot(h1_,h2_);

end