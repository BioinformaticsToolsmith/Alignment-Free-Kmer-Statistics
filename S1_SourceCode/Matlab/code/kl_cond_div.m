function div = kl_cond_div(h1, h2)
%   KL_COND_DIV computes the Kullback-Liebler Conditional divergence between the two
%   histograms

div1=0;
div2=0;
h1_sum=sum(h1);
h2_sum=sum(h2);

for i=1:4:size(h1,2)
    priorcount=sum(h1([i i+1 i+2 i+3]));
    priorcount2=sum(h2([i i+1 i+2 i+3]));
    innersum1=0;
    innersum2=0;
    for n=i:i+3
        % Direction 1
        h1cond=(h1(n)./priorcount);
        h2cond=(h2(n)./priorcount2);
        logterm=log(h1cond/h2cond);
        innersum1=innersum1+(h1cond*logterm);
        
        % Direction 2
        innersum2=innersum2+(h2cond*-logterm);
    end
    div1=div1+(priorcount/h1_sum)*innersum1;
    div2=div2+(priorcount2/h2_sum)*innersum2;
end
div = (div1+div2)/2;
end
