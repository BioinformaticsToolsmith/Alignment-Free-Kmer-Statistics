function newFeatures = feature_pairs(Features)
%   FEATURE_PAIRS computes the multiplicative statistical pairs
%
%
% [r, c] = size(Features);
% newFeatures=zeros(r, c * (c-1) / 2);
% 
% for n=1:c
%     f1=Features(:,n);
%     for m=n+1:c
%         f2=Features(:,m);
%         newFeatures(:,n) = f1.*f2;
%     end
% end
% 
% end

newFeatures=[];
stop=size(Features,2);
for n=1:stop
    f1=Features(:,n);
    for m=n+1:stop
        f2=Features(:,m);
        result=f1.*f2;
        newFeatures=[newFeatures result];
    end
end

end
