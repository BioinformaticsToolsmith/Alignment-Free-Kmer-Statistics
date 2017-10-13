function Ranks = review_paper_topNN(Table,Align_data,filename,Names,numFeature_end,k)
%   REVIEW_PAPER_TOPNN runs the Nearest Neighbor Experiment
%   -Table is the Statistics Table with each column representing one
%   statistic
%   -Align_data is the identity score table with sequence numbers
%   -filename is the name of the file where the results will be written
%   -Names is the Names of all the single/paired statistics in a cell
%   -numFeature_end is the number of the last statistic, use 33 to evaluate
%   single statistics only
%   -k is the K value for the K-nearest Neighbors, how many neighbors


start=1;
Features=[1:numFeature_end];
stop=size(Align_data,1);
num=stop-start+1;

Predict=Table;

numF=size(Features,2);
TopCorrect=zeros(numF,2);
Ranks=zeros(num,numF);

k_=k-1;
clear=stop-k_;

Ranks_align=[stop:-1:clear];
Features_good=[];

for i=1:numF
    
    feature=Predict(:,i);
    [~,I]=sort(feature);
    Ranks(:,i)=I;
    TopCorrect(i,1)=i;
    count=0;
    for j=stop:-1:clear
        if ismember(Ranks(j,i),Ranks_align)
            count=count+1;
        end
    end
    TopCorrect(i,2)=count;
    if count==size(Ranks_align,2)
        Features_good=[Features_good; i];
    end
end
[~,J]=sort(TopCorrect(:,2),'descend');
TopCorrect=TopCorrect(J,:);
writefile_names(filename,TopCorrect,1,Names,'KNN');

end