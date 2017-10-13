function d = pearson_coeff(h1, h2)
%   PEARSON_COEFF computes the Pearson Correlation Coefficient between the
%   two histograms
%
%

mu1=mean(h1);
mu2=mean(h2);

d_top=sum((h1-mu1).*(h2-mu2));
d_bottom= sqrt((sum((h1-mu1).^2)*sum((h2-mu2).^2)));
d=d_top/d_bottom;

end