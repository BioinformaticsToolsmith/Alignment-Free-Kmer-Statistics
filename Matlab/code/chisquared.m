function d = chisquared ( h1 , h2 )
% CHISQUARED computes the chi-squared distance for two histograms h1 and h2

d_item = ( h1 - h2 ).^2 ./  (h1 + h2 );
d_item(isnan(d_item) | isinf(d_item)) = 0;
d = sum(d_item);
end