function d = AFd(h1_2mer, h2_2mer, h1_1mer, h2_1mer)
%   AFD computes the Absolute Frequency Distance between two sequences
%   -h1_2mer and h2_2mer are the 2-mer histograms for each sequence
%   -h1_1mer and h2_1mer are the 1-mer histograms for each sequence
%

h1p=zeros(1,16);
h2p=zeros(1,16);

for i=1:16
    if i<=4
        h1p(i)=h1_2mer(i)/h1_1mer(1);
        h2p(i)=h2_2mer(i)/h2_1mer(1);
    elseif i<=8
        h1p(i)=h1_2mer(i)/h1_1mer(2);
        h2p(i)=h2_2mer(i)/h2_1mer(2);
    elseif i<=12
       h1p(i)=h1_2mer(i)/h1_1mer(3);
       h2p(i)=h2_2mer(i)/h2_1mer(3);
    else    
        h1p(i)=h1_2mer(i)/h1_1mer(4);
        h2p(i)=h2_2mer(i)/h2_1mer(4);
    end   
end
h_diff=h1p-h2p;

d= sum(((h_diff).*(1./((1+h_diff).^14))).^2);



end