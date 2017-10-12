function d = N2_r(h1, h2,lookup_r)
%   N2_R computes the N2_r similarity between the two histograms
%   -lookup_r is a lookup table that provides the indices of the reverse
%   k-mers

numKmers=size(h1,2);
k=log(numKmers) / log(4);

indices_r=zeros(numKmers,1);
for i=1:numKmers
    row=lookup_r(i,:);
    for j=1:k
       indices_r(i)=indices_r(i)+(row(j)*4^(k-j));
    end
end
indices_r=indices_r+1;
h1_=h1+h1(indices_r);
h2_=h2+h2(indices_r);
h1_z=(h1_-mean(h1_))/std(h1_);
h2_z=(h2_-mean(h2_))/std(h2_);
h1_zn=h1_z/norm(h1_z);
h2_zn=h2_z/norm(h2_z);
d=dot(h1_zn,h2_zn);

end