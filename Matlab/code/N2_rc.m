function d = N2_rc(h1, h2,lookup_rc)
%   N2_RC computes the N2_rc similarity between the two histograms
%   -lookup_rc is a lookup table that provides the indices of the reverse
%   complement
%   k-mers
numKmers=size(h1,2);
k=log(numKmers) / log(4);

indices_rc=zeros(numKmers,1);
for i=1:numKmers
    row=lookup_rc(i,:);
    for j=1:k
       indices_rc(i)=indices_rc(i)+(row(j)*4^(k-j));
    end
end
indices_rc=indices_rc+1;
h1_=h1+h1(indices_rc);
h2_=h2+h2(indices_rc);
h1_z=(h1_-mean(h1_))/std(h1_);
h2_z=(h2_-mean(h2_))/std(h2_);
h1_zn=h1_z/norm(h1_z);
h2_zn=h2_z/norm(h2_z);
d=dot(h1_zn,h2_zn);

end