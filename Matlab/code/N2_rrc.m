function d = N2_rrc(h1, h2,lookup_r,lookup_rc)
%   N2_RRC computes the N2_rrc similarity between the two histograms
%   -lookup_r is a lookup table that provides the indices of the reverse
%   k-mers
%   -lookup_rc is a lookup table that provides the indices of the reverse
%   complement k-mers
numKmers=size(h1,2);
k=log(numKmers) / log(4);

indices_r=zeros(numKmers,1);
indices_rc=zeros(numKmers,1);
for i=1:numKmers
    row=lookup_r(i,:);
    row2=lookup_rc(i,:);
    for j=1:k
       indices_r(i)=indices_r(i)+(row(j)*4^(k-j));
       indices_rc(i)=indices_rc(i)+(row2(j)*4^(k-j));
    end
end
indices_r=indices_r+1;
indices_rc=indices_rc+1;
h1_=h1+h1(indices_r)+h1(indices_rc);
h2_=h2+h2(indices_r)+h2(indices_rc);
h1_z=(h1_-mean(h1_))/std(h1_);
h2_z=(h2_-mean(h2_))/std(h2_);
h1_zn=h1_z/norm(h1_z);
h2_zn=h2_z/norm(h2_z);
d=dot(h1_zn,h2_zn);

end