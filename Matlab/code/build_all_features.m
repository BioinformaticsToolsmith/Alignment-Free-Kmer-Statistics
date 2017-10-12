function T= build_all_features(Align_table,Hist,Hist1,Hist2,start, stop,varargin)
%   BUILD_ALL_FEATURES builds all the statistics based on a given set 
%       of alignments, k-mer histograms, 1-mer histograms, and 2-mer histograms
%   Align_table(:,[2 3]) provides the sequence numbers for each statistic


    pack;
    Names= {'Hellinger' 'Manhattan' 'Euclidean' 'Chi-Squared' 'Normalized-Vectors' 'Harmonic-Mean' 'Jefferey-Div' 'K-Div' 'Pearson-Coeff' 'Squared-Chord' 'KL-Cond' 'Markov' 'Intersection' 'RRE-k-r' 'D2z' 'SimMM' 'EuclideanZ' 'EMD' 'Spearman' 'Jaccard' 'Lengthd' 'D2s' 'AFd' 'Mismatch' 'Canberra' 'Kulczynski1' 'Kulczynski2' 'Similarity-ratio' 'Jensen-Shannon' 'D2-star' 'N2r' 'N2rc' 'N2rrc'};
    f1=zeros(stop-start+1,1);
    f2=zeros(stop-start+1,1);
    f3=zeros(stop-start+1,1);
    f4=zeros(stop-start+1,1);
    f5=zeros(stop-start+1,1);
    f6=zeros(stop-start+1,1);
    f7=zeros(stop-start+1,1);
    f8=zeros(stop-start+1,1);
    f9=zeros(stop-start+1,1);
    f10=zeros(stop-start+1,1);
    f11=zeros(stop-start+1,1);
    f12=zeros(stop-start+1,1);
    f13=zeros(stop-start+1,1);
    f14=zeros(stop-start+1,1);
    f15=zeros(stop-start+1,1);
    f16=zeros(stop-start+1,1);
    f17=zeros(stop-start+1,1);
    f18=zeros(stop-start+1,1);
    f19=zeros(stop-start+1,1);
    f20=zeros(stop-start+1,1);
    f21=zeros(stop-start+1,1);
    f22=zeros(stop-start+1,1);
    f23=zeros(stop-start+1,1);
    f24=zeros(stop-start+1,1);
    f25=zeros(stop-start+1,1);
    f26=zeros(stop-start+1,1);
    f27=zeros(stop-start+1,1);
    f28=zeros(stop-start+1,1);
    f29=zeros(stop-start+1,1);
    f30=zeros(stop-start+1,1);
    f31=zeros(stop-start+1,1);
    f32=zeros(stop-start+1,1);
    f33=zeros(stop-start+1,1);
    
    numKmers=size(Hist,2);
    k=log(numKmers) / log(4);
    current=[1; 2; 3; 4];
    next=[];
    for j=2:k
        for h=1:size(current,1)
            for i=1:4
                next = [next; [current(h,:) i]];
            end
        end
    current = next;
    next = [];
    end
    lookup=current;
     lookup_=lookup-1;
            lookup_r=fliplr(lookup_);
            lookup_rc=lookup_r;

            for a=1:numKmers        %make reverse complement table
                for b=1:k
                    elt=lookup_rc(a,b);
                    switch(elt)
                        case 0
                            lookup_rc(a,b)=3;
                        case 1
                            lookup_rc(a,b)=2;
                        case 2
                            lookup_rc(a,b)=1;
                        case 3
                            lookup_rc(a,b)=0;
                        otherwise
                            error(message('Incorrect lookup table'));
                    end
                end
            end
    
    for n=start:stop
        h1=Hist(Align_table(n,2),:);
        h2=Hist(Align_table(n,3),:);
        h1_1mer=(Hist1(Align_table(n,2),:))-1;
        h2_1mer=(Hist1(Align_table(n,3),:))-1;
        h1_2mer=(Hist2(Align_table(n,2),:))-1;
        h2_2mer=(Hist2(Align_table(n,3),:))-1;
        h1_nonprob=h1-1;
        h2_nonprob=h2-1;
        f1(n,1)=hellinger(h1_nonprob,h2_nonprob);
        f2(n,1)=manhattan(h1_nonprob,h2_nonprob);
        f3(n,1)=euclidean(h1_nonprob,h2_nonprob);
        f4(n,1)=chisquared(h1_nonprob,h2_nonprob);
        f5(n,1)=normalized_vectors(h1_nonprob,h2_nonprob);
        f6(n,1)=harmonic_mean(h1,h2);
        f7(n,1)=jefferey_divergence(h1,h2);
        f8(n,1)=k_divergence(h1,h2);
        f9(n,1)=pearson_coeff(h1_nonprob,h2_nonprob);
        f10(n,1)=squaredchord(h1_nonprob,h2_nonprob);
        f11(n,1)=kl_cond_div(h1,h2);
        model1= markov(h1,h2);
        model2= markov(h2, h1);
        f12(n,1)=(model1+model2)/2;
        f13(n,1)=intersection(h1_nonprob,h2_nonprob);
        f14(n,1)=rre_k_r(h1,h2);
        f15(n,1)=D2z(h1_nonprob,h2_nonprob);
        f16(n,1)=SimMM(h1,h2);
        f17(n,1)=euclideanZ(h1_nonprob,h2_nonprob);                %standardized euclidean
        %f18(n,1)=pdist_(h1,h2,'emd');        %changing EMD
        f18(n,1)=EMD(h1,h2);
        f19(n,1)=pdist([h1; h2],'spearman');
        f20(n,1)=pdist([h1; h2],'jaccard');
        f21(n,1)=lengthd(h1_nonprob,h2_nonprob);
        f22(n,1)=D2s(h1_nonprob,h2_nonprob,h1_1mer,h2_1mer,lookup);
        f23(n,1)=AFd(h1_2mer,h2_2mer,h1_1mer,h2_1mer);
        f24(n,1)=mismatch(h1_nonprob,h2_nonprob);
        f25(n,1)=canberra(h1,h2);
        f26(n,1)=kulczynski1(h1,h2);
        f27(n,1)=kulczynski2(h1_nonprob,h2_nonprob);
        f28(n,1)=similarity_ratio(h1_nonprob,h2_nonprob);
        f29(n,1)=jensen_shannon(h1,h2);
        f30(n,1)=D2_star(h1_nonprob,h2_nonprob,h1_1mer,h2_1mer,lookup);
        f31(n,1)=N2_r(h1_nonprob,h2_nonprob,lookup_r);
        f32(n,1)=N2_rc(h1_nonprob,h2_nonprob,lookup_rc);
        f33(n,1)=N2_rrc(h1_nonprob,h2_nonprob,lookup_r,lookup_rc);
        n
    end
    
    T=[f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 f31 f32 f33];
    
    fix=[1 2 3 4 7 8 10 11 14 17 18 19 20 21 23 24 25 26 29];

    numF=size(T,2);
    for b=1:numF
        T(:,b)=normalize(T(:,b));
    end
    
    for i=1:numF
        if ismember(i,fix)
            T(:,i)=1-T(:,i);
        end   
    end
    
    if nargin>6
        Exp_name=[' - ' varargin];
    figure();
    for j=1:numF
        subplot(6,6,j);
        plot(Align_table(:,1),T(:,j),'.');
        title([Names{j} Exp_name] );
    end
    
    end
    
    T=[T T.^2];
    T=[T feature_pairs(T)];
    
end