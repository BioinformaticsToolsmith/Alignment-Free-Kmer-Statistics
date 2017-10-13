function probs= build_cond(hist)
%BUILD_COND builds a conditional probability vector for a given histogram

numkmers=size(hist,2);
probs=zeros(1,numkmers);
    for i=1:4:numkmers
        sumOfFour = sum(hist([i i+1 i+2 i+3]));
        for k=i+0:i+3
            probs(1,k)=(hist(k)/sumOfFour);
        end
    end
end

