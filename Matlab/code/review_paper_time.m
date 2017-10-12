
function Times= review_paper_time(filename,Align_table,Hist,Hist1, Hist2,start,stop)
%   REVIEW_PAPER_TIME times how long the statistics take to evaluate each
%   of the sequence pairs

Times=zeros(33,2);
Times(:,1)=[1:33];
Names={'Hellinger' 'Manhattan' 'Euclidean' 'Chi-Squared' 'Norm-Vectors' 'Harmonic-Mean' 'Jefferey-Div' 'K-Div' 'Pearson-Coeff' 'Squared-Chord' 'KL-Cond' 'Markov' 'Intersection' 'RRE-k-r' 'D2z' 'SimMM' 'EuclideanZ' 'EMD' 'Spearman' 'Jaccard' 'Lengthd' 'D2s' 'AFd' 'Mismatch' 'Canberra' 'Kulczynski1' 'Kulczynski2' 'Similarity-ratio' 'Jensen-Shannon' 'D2-star' 'N2r' 'N2rc' 'N2rrc'};

Predict=zeros(stop,1);
for f=1:33
    
    switch(f)
        case 1
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=hellinger(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 2
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=manhattan(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 3
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=euclidean(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 4
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=chisquared(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 5
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=normalized_vectors(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 6
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=harmonic_mean(h1,h2);
            end
            toc;
            Times(f,2)=toc;
        case 7
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=jefferey_divergence(h1,h2);
            end
            toc;
            Times(f,2)=toc;
        case 8
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=k_divergence(h1,h2);
            end
            toc;
            Times(f,2)=toc;
        case 9
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=pearson_coeff(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 10
           
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=squaredchord(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 11
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=kl_cond_div(h1,h2);
            end
            toc;
            Times(f,2)=toc;
        case 12
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                model1= markov(h1,h2);
                model2= markov(h2, h1);
                Predict(n)=(model1+model2)/2;
            end
            toc;
            Times(f,2)=toc;
        case 13
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=intersection(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 14
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=rre_k_r(h1,h2);
            end
            toc;
            Times(f,2)=toc;
        case 15
           
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=D2z(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 16
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=SimMM(h1,h2);
            end
            toc;
            Times(f,2)=toc;
        case 17
           
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=euclideanZ(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 18
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=EMD(h1,h2);
            end
            toc;
            Times(f,2)=toc;
        case 19
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=pdist([h1; h2],'spearman');
            end
            toc;
            Times(f,2)=toc;
        case 20
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=pdist([h1; h2],'jaccard');
            end
            toc;
            Times(f,2)=toc;
        case 21
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=lengthd(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 22
            tic;
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
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=D2s(h1_nonprob,h2_nonprob,h1_1mer,h2_1mer,lookup);
            end
            toc;
            Times(f,2)=toc;
        case 23
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=AFd(h1_2mer,h2_2mer,h1_1mer,h2_1mer);
            end
            toc;
            Times(f,2)=toc;
        case 24
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=mismatch(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 25
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=canberra(h1,h2);
            end
            toc;
            Times(f,2)=toc;
        case 26
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=kulczynski1(h1,h2);
            end
            toc;
            Times(f,2)=toc;
        case 27
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=kulczynski2(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        case 28
            
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=similarity_ratio(h1_nonprob,h2_nonprob);
            end
            toc;
            Times(f,2)=toc;
        
        case 29    
            tic;
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=jensen_shannon(h1,h2);
            end
            toc;
            Times(f,2)=toc;
        case 30    
            tic;
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
            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=D2_star(h1_nonprob,h2_nonprob,h1_1mer,h2_1mer,lookup);
            end
            toc;
            Times(f,2)=toc;
        case 31   
            tic;
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

            for n=start:stop
                h1=Hist(Align_table(n,2),:);
                h2=Hist(Align_table(n,3),:);
                h1_1mer=(Hist1(Align_table(n,2),:))-1;
                h2_1mer=(Hist1(Align_table(n,3),:))-1;
                h1_2mer=(Hist2(Align_table(n,2),:))-1;
                h2_2mer=(Hist2(Align_table(n,3),:))-1;
                h1_nonprob=h1-1;
                h2_nonprob=h2-1;
                Predict(n)=N2_r(h1_nonprob,h2_nonprob,lookup_r);
            end
            toc;
            Times(f,2)=toc;
        case 32   
            tic;
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
            lookup_rc=fliplr(lookup_);
            for a=1:numKmers
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
                Predict(n)=N2_rc(h1_nonprob,h2_nonprob,lookup_rc);
            end
            toc;
            Times(f,2)=toc;
        case 33  
            tic;
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
                Predict(n)=N2_rrc(h1_nonprob,h2_nonprob,lookup_r,lookup_rc);
            end
            toc;
            Times(f,2)=toc;
        otherwise
            error(message('incorrect number of features'));
    end
end
[~,I]=sort(Times(:,2),'descend');
Times=Times(I,:);

figure();
F=barh(Times(:,2),.4);
set(gca,'ytick',1:33,'yticklabel',Names(Times(:,1)))
xlabel('Time (seconds)');
xlim([0 (Times(1,2)+3)]);
ylabel('Statistics');
ylim([0 34]);
xt=Times(:,2);
yt=[1:33];
ytxt=num2str(Times(:,2),'%.3f');
text(xt,yt,ytxt,'fontweight','bold');
set(gca, 'FontSize', 12);
writefile_names(filename,Times,1,Names,'Times');

end