function d = D2s(h1, h2, h1_1mer, h2_1mer,lookup)
% D2S computes the D2s distance for two histograms h1 and h2 as
% well as their corresponding 1-mer histograms
% lookup is a table provided for k-mer indexing
numKmers=size(h1,2);
k=log(numKmers) / log(4);


% Converting the mono-mer table to probabilities distributions
count1=sum(h1_1mer);
count2=sum(h2_1mer);
p1=[h1_1mer(1)/count1 h1_1mer(2)/count1 h1_1mer(3)/count1 h1_1mer(4)/count1];
p2=[h2_1mer(1)/count2 h2_1mer(2)/count2 h2_1mer(3)/count2 h2_1mer(4)/count2];

% 
h1_=zeros(1,numKmers);
h2_=zeros(1,numKmers);

for i=1:numKmers
    p_1i=1;
    p_2i=1;
    indices=lookup(i,:);
    for j=1:k
        p_1i=p_1i*p1(indices(j));
        p_2i=p_2i*p2(indices(j));
    end
    h1_(i)=h1(i)-(count1)*p_1i;
    h2_(i)=h2(i)-(count2)*p_2i;
end

d= sum((h1_.*h2_) ./ sqrt((h1_).^2 + (h2_).^2));


end