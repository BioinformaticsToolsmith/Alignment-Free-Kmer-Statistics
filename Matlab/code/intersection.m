function d = intersection(h1, h2)
%   INTERSECTION computes the intersection similarity between the two
%   histograms


d=(2*sum(min(h1,h2)))/(sum(h1)+sum(h2));

end