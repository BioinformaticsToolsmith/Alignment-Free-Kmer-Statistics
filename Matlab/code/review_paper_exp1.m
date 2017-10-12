
function SensSpec = review_paper_exp1(filename,Data,Hist,Hist1,Hist2,min_lower,min_upper,threshold_lower,threshold_upper,sample_bad,sample_good)
%   REVIEW_PAPER_EXP1 runs the Sensitivity/Specificity Experiment
%   min_lower to threshold_lower specifies the range of the identity
%       scores where the negative set will be randomly selected
%   min_upper to threshold_upper specifies the range of the identity
%       scores where the positive set will be randomly selected

Features=[1:33];

Names={'Hellinger' 'Manhattan' 'Euclidean' 'Chi-Squared' 'Normalized-Vectors' 'Harmonic-Mean' 'Jefferey-Div' 'K-Div' 'Pearson-Coeff' 'Squared-Chord' 'KL-Cond' 'Markov' 'Intersection' 'RRE-k-r' 'D2z' 'SimMM' 'EuclideanZ' 'EMD' 'Spearman' 'Jaccard' 'Lengthd' 'D2s' 'AFd' 'Mismatch' 'Canberra' 'Kulczynski1' 'Kulczynski2' 'Similarity-ratio' 'Jensen-Shannon' 'D2-star' 'N2r' 'N2rc' 'N2rrc'};
rng('default');
%Data = Data(Data(:,2) == 1, :); 
Align_bad_set1=Data(Data(:,1)>=min_lower,:);
Align_bad_set2=Data(Data(:,1)<=min_upper,:);
Align_bad_set=intersect(Align_bad_set1, Align_bad_set2,'rows');
size1=size(Align_bad_set,1);
disp(['review_paper_exp1 line 18 - size 1 is: ' num2str(size1)]);
if sample_bad>size1
    error(message('not enough data points in negative set'));
else
    y_bad=randsample(size1,sample_bad);
end
Align_bad=Align_bad_set(y_bad,:);

Align_good_set1=Data(Data(:,1)>=threshold_lower,:);
Align_good_set2=Data(Data(:,1)<=threshold_upper,:);
Align_good_set=intersect(Align_good_set1, Align_good_set2,'rows');
size2=size(Align_good_set,1);
disp(['review_paper_exp1 line 30 - size 2 is: ' num2str(size2)]);

if sample_good>size2
    error(message('not enough data points in positive set'));
else
    y_good=randsample(size2,sample_good);
end
Align_good=Align_good_set(y_good,:);
Align=[Align_bad; Align_good];

num=sample_bad+sample_good;

if num ~= size(Align, 1)
    error(message('sample_bad+sample_good != size(Align,1)'));
end

Table=build_all_features(Align,Hist,Hist1,Hist2,1,num);
Bad_seqs=Align([1:sample_bad],[2:3]);
Good_seqs=Align([(sample_bad+1):num],[2:3]);
Seqs=[Bad_seqs; Good_seqs];
SensSpec=zeros(size(Features,2),2);
SensSpec(:,1)=Features;
disp(['Balanced the data']);

for i=Features
    disp(['Feature ' int2str(i) ' out of ' int2str(size(Features, 2))]);
    correctgood=0;
    correctbad=0;
    incorrect=[];
    feature=[Table(:,i) Seqs];
    [B, I] =sort(feature(:,1));
    feature=feature(I,:);
    % Commented on 10/4 by Hani and Ben
    % Align_compare=zeros(num,1);
    for j=1:num
        seqs=feature(j,[2:3]);
        row=find(ismember(Align(:,[2 3]),seqs,'rows'));
        if j<=sample_bad
            if Align(row,1)<threshold_lower
                correctbad=correctbad+1;
                
            else
                incorrect=[incorrect j];
            end
        else
            if Align(row,1)>=threshold_lower
               correctgood=correctgood+1;
            else
                incorrect=[incorrect j];
            end
        end
        % Commented on 10/4 by Hani and Ben
        % Align_compare(j)=Align(row,1);
        
    end
    
    SensSpec(i,2)=correctgood/sample_good;
    SensSpec(i,3)=correctbad/sample_bad;
    
end
[B,J]=sort(SensSpec(:,2),'descend');
SensSpec=SensSpec(J,:);

F=figure();
hold on
for i = 1:length(SensSpec(:,2))
    h=bar(i,SensSpec(i,2),.7);
    if SensSpec(i,2) < 0.9
        set(h,'FaceColor','k');
    elseif SensSpec(i,2) <= 0.99
        set(h,'FaceColor','b');
    else
        set(h,'FaceColor','r');
    end
end
hold off

set(gca,'xtick',1:33,'xticklabel',Names(SensSpec(:,1)),'XTickLabelRotation',90);
ylabel('Sensitivity');
set(gca, 'FontSize', 14);
saveas(F, [filename '/Sensitivity'], 'epsc');    

F2=figure();
hold on
for i = 1:length(SensSpec(:,3))
    h=bar(i,SensSpec(i,3),.7);
    if SensSpec(i,3) < 0.9
        set(h,'FaceColor','k');
    elseif SensSpec(i,3) < 1
        set(h,'FaceColor','b');
    else
        set(h,'FaceColor','r');
    end
end
hold off

set(gca,'xtick',1:33,'xticklabel',Names(SensSpec(:,1)),'XTickLabelRotation',90);
ylabel('Specificity');
set(gca, 'FontSize', 14);
saveas(F2, [filename '/Specificity'], 'epsc');  

Y_p27=discretize(Data,[0:.05:1]);
T_p27=makebins(Data,Y_p27,20);
sizes_p27=zeros(20,1);
for n=1:20
sizes_p27(n)=size(T_p27{n},1);
end

 F=bar(sizes_p27);
 xlabel('Alignment Identity Score Bins (each 5%)');
 ylabel('Number of Sequence Pairs');
 set(gca,'xtick',1:20,'xticklabel',{'0-5' '5-10' '10-15' '15-20' '20-25' '25-30' '30-35' '35-40' '40-45' '45-50' '50-55' '55-60' '60-65' '65-70' '70-75' '75-80' '80-85' '85-90' '90-95' '95-100'},'XTickLabelRotation',90);
 xlim([0 21]);
 saveas(F, [filename '/Distribution'], 'epsc'); 


end