function Review_Paper_Results(output_dir)
%   'Review_Paper_Results' duplicates the
%   results from the Review paper
%   output_dir is a string of the path to the output_dir
%   -this is the driver program for running each of the experiments with
%   the paramaters used in the review
%   -to run the evaluation benchmark, see EVALUATE and REVIEW_PAPER

warning ('off','all');

load('../Data/Example_Data.mat');
directory=output_dir;
mkdir(directory);
mkdir ([directory '/exp1']);
mkdir ([directory '/exp2']);
mkdir ([directory '/exp2/Fig']);
mkdir ([directory '/exp3']);
mkdir ([directory '/exp3/p27']);
mkdir ([directory '/exp3/microbial']);
mkdir ([directory '/exp3/microbial/kmersize3']);
mkdir ([directory '/exp3/microbial/kmersize2']);

disp('Running Sensitivity Experiment');
review_paper_exp1([directory '/exp1'],Align_p27_1seq,Hist6_p27,Hist1_p27,Hist2_p27,0,.7,.7,1,109,26);

disp('Building Table for Linear Correlation experiment');
numData=size(AlignNGS_data_plotting,1);
Table=build_all_features(AlignNGS_data_plotting,HistNGS,HistNGS_1,HistNGS_2,1,numData,'R^2,Micrbiome');

numFeatures=size(Table,2);
Names_singles={'Hellinger' 'Manhattan' 'Euclidean' 'Chi-Squared' 'Normalized-Vectors' 'Harmonic-Mean' 'Jefferey-Div' 'K-Div' 'Pearson-Coeff' 'Squared-Chord' 'KL-Cond' 'Markov' 'Intersection' 'RRE-k-r' 'D2z' 'SimMM' 'EuclideanZ' 'EMD' 'Spearman' 'Jaccard' 'Lengthd' 'D2s' 'AFd' 'Mismatch' 'Canberra' 'Kulczynski1' 'Kulczynski2' 'Similarity-ratio' 'Jensen-Shannon' 'D2-star' 'N2r' 'N2rc' 'N2rrc'};
Names=make_feature_names(Names_singles);
cut60=numData;
cut70=numData;
cut80=numData;
cut90=numData;
for i=1:numData
    current=AlignNGS_data_plotting(i,1);
    if current>.9 && current<AlignNGS_data_plotting(cut90,1)
        cut90=i;
        continue
    end
    if current>.8 && current<AlignNGS_data_plotting(cut80,1)
        cut80=i;
        continue
    end
    if current>.7 && current<AlignNGS_data_plotting(cut70,1)
        cut70=i;
        continue
    end
    if current>.6 && current<AlignNGS_data_plotting(cut60,1)
        cut60=i;
        continue
    end
end

disp('Linear Correlation: 0% threshold:');
review_paper_plot([directory '/exp2'],[directory '/exp2/Rsquared_all'],Table,AlignNGS_data_plotting,Names,1,0,numFeatures);
disp('Linear Correlation: 60% threshold:');
review_paper_plot([directory '/exp2'],[directory '/exp2/Rsquared_60'],Table,AlignNGS_data_plotting,Names,cut60,.6,numFeatures);
disp('Linear Correlation: 70% threshold:');
review_paper_plot([directory '/exp2'],[directory '/exp2/Rsquared_70'],Table,AlignNGS_data_plotting,Names,cut70,.7,numFeatures);
disp('Linear Correlation: 80% threshold:');
review_paper_plot([directory '/exp2'],[directory '/exp2/Rsquared_80'],Table,AlignNGS_data_plotting,Names,cut80,.8,numFeatures);
disp('Linear Correlation: 90% threshold:');
review_paper_plot([directory '/exp2'],[directory '/exp2/Rsquared_90'],Table,AlignNGS_data_plotting,Names,cut90,.9,numFeatures);

disp('Building p27 Table for KNN experiment');
Table1=build_all_features(Align_p27_1seq,Hist6_p27,Hist1_p27,Hist2_p27,1,size(Align_p27_1seq,1),'SensSpec,p27');
disp('K=1, singles:');
review_paper_topNN(Table1,Align_p27_1seq,[directory '/exp3/p27/K=1,s'],Names,33,1);
disp('K=1, pairs:');
review_paper_topNN(Table1,Align_p27_1seq,[directory '/exp3/p27/K=1,p'],Names,2211,1);
disp('K=5, singles:');
review_paper_topNN(Table1,Align_p27_1seq,[directory '/exp3/p27/K=5,s'],Names,33,5);
disp('K=5, pairs:');
review_paper_topNN(Table1,Align_p27_1seq,[directory '/exp3/p27/K=5,p'],Names,2211,5);
disp('K=10, singles:');
review_paper_topNN(Table1,Align_p27_1seq,[directory '/exp3/p27/K=10,s'],Names,33,10);
disp('K=10, pairs:');
review_paper_topNN(Table1,Align_p27_1seq,[directory '/exp3/p27/K=10,p'],Names,2211,10);

disp('Building Microbiome Table for KNN experiment');
Table2=build_all_features(AlignNGS_KNN,HistNGS,HistNGS_1,HistNGS_2,1,size(AlignNGS_KNN,1),'KNN,Microbiome');

disp('K=1, singles:');
review_paper_topNN(Table2,AlignNGS_KNN,[directory '/exp3/microbial/K=1,s'],Names,33,1);
disp('K=1, pairs:');
review_paper_topNN(Table2,AlignNGS_KNN,[directory '/exp3/microbial/K=1,p'],Names,2211,1);
disp('K=5, singles:');
review_paper_topNN(Table2,AlignNGS_KNN,[directory '/exp3/microbial/K=5,s'],Names,33,5);
disp('K=5, pairs:');
review_paper_topNN(Table2,AlignNGS_KNN,[directory '/exp3/microbial/K=5,p'],Names,2211,5);
disp('K=10, singles:');
review_paper_topNN(Table2,AlignNGS_KNN,[directory '/exp3/microbial/K=10,s'],Names,33,10);
disp('K=10, pairs:');
review_paper_topNN(Table2,AlignNGS_KNN,[directory '/exp3/microbial/K=10,p'],Names,2211,10);

disp('Building Table for KNN experiment with k-mer size=3');
Table3=build_all_features(AlignNGS_KNN,HistNGS_3,HistNGS_1,HistNGS_2,1,size(AlignNGS_KNN,1),'KNN,k minus 1');

disp('K=1, singles:');
review_paper_topNN(Table3,AlignNGS_KNN,[directory '/exp3/microbial/kmersize3/K=1,s'],Names,33,1);
disp('K=1, pairs:');
review_paper_topNN(Table3,AlignNGS_KNN,[directory '/exp3/microbial/kmersize3/K=1,p'],Names,2211,1);
disp('K=5, singles:');
review_paper_topNN(Table3,AlignNGS_KNN,[directory '/exp3/microbial/kmersize3/K=5,s'],Names,33,5);
disp('K=5, pairs:');
review_paper_topNN(Table3,AlignNGS_KNN,[directory '/exp3/microbial/kmersize3/K=5,p'],Names,2211,5);
disp('K=10, singles:');
review_paper_topNN(Table3,AlignNGS_KNN,[directory '/exp3/microbial/kmersize3/K=10,s'],Names,33,10);
disp('K=10, pairs:');
review_paper_topNN(Table3,AlignNGS_KNN,[directory '/exp3/microbial/kmersize3/K=10,p'],Names,2211,10);

disp('Building Table for KNN experiment with k-mer size=2');
Table4=build_all_features(AlignNGS_KNN,HistNGS_2,HistNGS_1,HistNGS_2,1,size(AlignNGS_KNN,1),'KNN,k minus 2');

disp('K=1, singles:');
review_paper_topNN(Table4,AlignNGS_KNN,[directory '/exp3/microbial/kmersize2/K=1,s'],Names,33,1);
disp('K=1, pairs:');
review_paper_topNN(Table4,AlignNGS_KNN,[directory '/exp3/microbial/kmersize2/K=1,p'],Names,2211,1);
disp('K=5, singles:');
review_paper_topNN(Table4,AlignNGS_KNN,[directory '/exp3/microbial/kmersize2/K=5,s'],Names,33,5);
disp('K=5, pairs:');
review_paper_topNN(Table4,AlignNGS_KNN,[directory '/exp3/microbial/kmersize2/K=5,p'],Names,2211,5);
disp('K=10, singles:');
review_paper_topNN(Table4,AlignNGS_KNN,[directory '/exp3/microbial/kmersize2/K=10,s'],Names,33,10);
disp('K=10, pairs:');
review_paper_topNN(Table4,AlignNGS_KNN,[directory '/exp3/microbial/kmersize2/K=10,p'],Names,2211,10);

disp('Time test:');
review_paper_time([directory '/Times'],AlignNGS_data_official,HistNGS,HistNGS_1,HistNGS_2,1,size(AlignNGS_data_official,1));

end