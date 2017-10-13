
function REVIEW_PAPER(Align_data,Histk,Histk_,Histk__,Histkp,Histkpp,Hist1,Hist2,path)

%REVIEW_PAPER runs all experiments, analyzes each of the statistics by using histograms
%on the alignment identity score file
% m=number of sequences
% n=number of sequence comparisons
%-Align_data is an n x 3 matrix where the first column is the alignment
%    identity score and columns 2 and 3 are the sequence numbers being compared
%    Ex:[0.8597 1 26]; sequence 1 and sequence 26 in the Histogram file have an
%    alignment identity score of 0.8597
%
%-Histk is an m x 4^k matrix where each row represents the k-mer counts for
%    that particular sequence
%
%-Hist1 is Histk when k=1, each row is only 4^1 k-mer count values
%
%-Hist2 is Histk when k=2, each row is only 4^2 k-mer count values
%
%-path is the output path that locates the target directory
%
%-results can be found in the target directory

directory=[path '/Matlab_EVAL'];
mkdir(directory);
mkdir ([directory '/exp1']);
mkdir ([directory '/exp2']);
mkdir ([directory '/exp2/Fig']);
mkdir ([directory '/exp3']);
mkdir ([directory '/exp3/kmersize_minus_1']);
mkdir ([directory '/exp3/kmersize_minus_2']);
mkdir ([directory '/exp3/kmersize_plus_1']);
mkdir ([directory '/exp3/kmersize_plus_2']);


numData=size(Align_data,1);
disp(numData);
[~,I]=sort(Align_data(:,1));
Align_data=Align_data(I,:);
I = find(Align_data(:, 2) == 1 & Align_data(:,1) ~= 1);
% Reduced_data=Align_data(Align_data(:, 2) == 1 & Align_data(:,1) ~= 1,:);
Reduced_data=Align_data(I,:);

Align_data = Reduced_data;
numData=size(Align_data, 1);

disp('Building Statisitcs ...');
Table=build_all_features(Align_data,Histk,Hist1,Hist2,1,numData);
numFeatures=size(Table,2);
Names_singles={'Hellinger' 'Manhattan' 'Euclidean' 'Chi-Squared' 'Normalized-Vectors' 'Harmonic-Mean' 'Jefferey-Div' 'K-Div' 'Pearson-Coeff' 'Squared-Chord' 'KL-Cond' 'Markov' 'Intersection' 'RRE-k-r' 'D2z' 'SimMM' 'EuclideanZ' 'EMD' 'Spearman' 'Jaccard' 'Lengthd' 'D2s' 'AFd' 'Mismatch' 'Canberra' 'Kulczynski1' 'Kulczynski2' 'Similarity-ratio' 'Jensen-Shannon' 'D2-star' 'N2r' 'N2rc' 'N2rrc'};
Names=make_feature_names(Names_singles);     %creates paired statistic names from the singles
cut60=numData;
cut70=numData;
cut80=numData;
cut90=numData;
for i=1:numData
    current=Align_data(i,1);
    if current>.9 && current<Align_data(cut90,1)
        cut90=i;
        continue
    end
    if current>.8 && current<Align_data(cut80,1)
        cut80=i;
        continue
    end
    if current>.7 && current<Align_data(cut70,1)
        cut70=i;
        continue
    end
    if current>.6 && current<Align_data(cut60,1)
        cut60=i;
        continue
    end
end

disp([cut60 cut70 cut80 cut90]);
%good_sample=min(numData-cut80,cut80-1);
cutoff = 0.9;
good_sample=size(find(Reduced_data(:,1) >= cutoff), 1);
bad_sample=good_sample;
disp([bad_sample good_sample]);
review_paper_exp1([directory '/exp1'],Reduced_data,Histk,Hist1,Hist2,0,cutoff,cutoff,1,bad_sample,good_sample);

disp('Linear Correlation: 0% threshold ...');
review_paper_plot([directory '/exp2'],[directory '/exp2/Rsquared_all'],Table,Align_data,Names,1,0,numFeatures);
disp('Linear Correlation: 60% threshold ...');
review_paper_plot([directory '/exp2'],[directory '/exp2/Rsquared_60'],Table,Align_data,Names,cut60,.6,numFeatures);
disp('Linear Correlation: 70% threshold ...');
review_paper_plot([directory '/exp2'],[directory '/exp2/Rsquared_70'],Table,Align_data,Names,cut70,.7,numFeatures);
disp('Linear Correlation: 80% threshold ...');
review_paper_plot([directory '/exp2'],[directory '/exp2/Rsquared_80'],Table,Align_data,Names,cut80,.8,numFeatures);
disp('Linear Correlation: 90% threshold ...');
review_paper_plot([directory '/exp2'],[directory '/exp2/Rsquared_90'],Table,Align_data,Names,cut90,.9,numFeatures);

disp('Running KNN experiment ...')
disp('K=1, singles:');
review_paper_topNN(Table,Reduced_data,[directory '/exp3/K=1,s'],Names,33,1);
disp('K=1, pairs:');
review_paper_topNN(Table,Reduced_data,[directory '/exp3/K=1,p'],Names,2211,1);
disp('K=5, singles:');
review_paper_topNN(Table,Reduced_data,[directory '/exp3/K=5,s'],Names,33,5);
disp('K=5, pairs:');
review_paper_topNN(Table,Reduced_data,[directory '/exp3/K=5,p'],Names,2211,5);
disp('K=10, singles:');
review_paper_topNN(Table,Reduced_data,[directory '/exp3/K=10,s'],Names,33,10);
disp('K=10, pairs:');
review_paper_topNN(Table,Reduced_data,[directory '/exp3/K=10,p'],Names,2211,10);

% disp('Running KNN experiment ...')
% disp('K=1, singles:');
% review_paper_topNN(Table(I,:),Reduced_data,[directory '/exp3/K=1,s'],Names,33,1);
% disp('K=1, pairs:');
% review_paper_topNN(Table(I,:),Reduced_data,[directory '/exp3/K=1,p'],Names,2211,1);
% disp('K=5, singles:');
% review_paper_topNN(Table(I,:),Reduced_data,[directory '/exp3/K=5,s'],Names,33,5);
% disp('K=5, pairs:');
% review_paper_topNN(Table(I,:),Reduced_data,[directory '/exp3/K=5,p'],Names,2211,5);
% disp('K=10, singles:');
% review_paper_topNN(Table(I,:),Reduced_data,[directory '/exp3/K=10,s'],Names,33,10);
% disp('K=10, pairs:');
% review_paper_topNN(Table(I,:),Reduced_data,[directory '/exp3/K=10,p'],Names,2211,10);


if(~isempty(Histk_))
disp('Building Table for KNN experiment with k-mer size reduced by 1 ...');
Table3=build_all_features(Reduced_data,Histk_,Hist1,Hist2,1,size(Reduced_data,1));

disp('K=1, singles:');
review_paper_topNN(Table3,Reduced_data,[directory '/exp3/kmersize_minus_1/K=1,s'],Names,33,1);
disp('K=1, pairs:');
review_paper_topNN(Table3,Reduced_data,[directory '/exp3/kmersize_minus_1/K=1,p'],Names,2211,1);
disp('K=5, singles:');
review_paper_topNN(Table3,Reduced_data,[directory '/exp3/kmersize_minus_1/K=5,s'],Names,33,5);
disp('K=5, pairs:');
review_paper_topNN(Table3,Reduced_data,[directory '/exp3/kmersize_minus_1/K=5,p'],Names,2211,5);
disp('K=10, singles:');
review_paper_topNN(Table3,Reduced_data,[directory '/exp3/kmersize_minus_1/K=10,s'],Names,33,10);
disp('K=10, pairs:');
review_paper_topNN(Table3,Reduced_data,[directory '/exp3/kmersize_minus_1/K=10,p'],Names,2211,10);

end

if(~isempty(Histk__))
disp('Building Table for KNN experiment with k-mer size reduced by 2 ...');
Table4=build_all_features(Reduced_data,Histk__,Hist1,Hist2,1,size(Reduced_data,1));

disp('K=1, singles:');
review_paper_topNN(Table4,Reduced_data,[directory '/exp3/kmersize_minus_2/K=1,s'],Names,33,1);
disp('K=1, pairs:');
review_paper_topNN(Table4,Reduced_data,[directory '/exp3/kmersize_minus_2/K=1,p'],Names,2211,1);
disp('K=5, singles:');
review_paper_topNN(Table4,Reduced_data,[directory '/exp3/kmersize_minus_2/K=5,s'],Names,33,5);
disp('K=5, pairs:');
review_paper_topNN(Table4,Reduced_data,[directory '/exp3/kmersize_minus_2/K=5,p'],Names,2211,5);
disp('K=10, singles:');
review_paper_topNN(Table4,Reduced_data,[directory '/exp3/kmersize_minus_2/K=10,s'],Names,33,10);
disp('K=10, pairs:');
review_paper_topNN(Table4,Reduced_data,[directory '/exp3/kmersize_minus_2/K=10,p'],Names,2211,10);

end
if(~isempty(Histkp))
disp('Building Table for KNN experiment with k-mer size increased by 1 ...');
Table5=build_all_features(Reduced_data,Histkp,Hist1,Hist2,1,size(Reduced_data,1));

disp('K=1, singles:');
review_paper_topNN(Table5,Reduced_data,[directory '/exp3/kmersize_plus_1/K=1,s'],Names,33,1);
disp('K=1, pairs:');
review_paper_topNN(Table5,Reduced_data,[directory '/exp3/kmersize_plus_1/K=1,p'],Names,2211,1);
disp('K=5, singles:');
review_paper_topNN(Table5,Reduced_data,[directory '/exp3/kmersize_plus_1/K=5,s'],Names,33,5);
disp('K=5, pairs:');
review_paper_topNN(Table5,Reduced_data,[directory '/exp3/kmersize_plus_1/K=5,p'],Names,2211,5);
disp('K=10, singles:');
review_paper_topNN(Table5,Reduced_data,[directory '/exp3/kmersize_plus_1/K=10,s'],Names,33,10);
disp('K=10, pairs:');
review_paper_topNN(Table5,Reduced_data,[directory '/exp3/kmersize_plus_1/K=10,p'],Names,2211,10);

end
if(~isempty(Histkpp))
disp('Building Table for KNN experiment with k-mer size increased by 2 ...');
Table6=build_all_features(Reduced_data,Histkpp,Hist1,Hist2,1,size(Reduced_data,1));

disp('K=1, singles:');
review_paper_topNN(Table6,Reduced_data,[directory '/exp3/kmersize_plus_2/K=1,s'],Names,33,1);
disp('K=1, pairs:');
review_paper_topNN(Table6,Reduced_data,[directory '/exp3/kmersize_plus_2/K=1,p'],Names,2211,1);
disp('K=5, singles:');
review_paper_topNN(Table6,Reduced_data,[directory '/exp3/kmersize_plus_2/K=5,s'],Names,33,5);
disp('K=5, pairs:');
review_paper_topNN(Table6,Reduced_data,[directory '/exp3/kmersize_plus_2/K=5,p'],Names,2211,5);
disp('K=10, singles:');
review_paper_topNN(Table6,Reduced_data,[directory '/exp3/kmersize_plus_2/K=10,s'],Names,33,10);
disp('K=10, pairs:');
review_paper_topNN(Table6,Reduced_data,[directory '/exp3/kmersize_plus_2/K=10,p'],Names,2211,10);

end
disp('Time test:');
review_paper_time([directory '/Times'],Align_data,Histk,Hist1,Hist2,1,numData);

end