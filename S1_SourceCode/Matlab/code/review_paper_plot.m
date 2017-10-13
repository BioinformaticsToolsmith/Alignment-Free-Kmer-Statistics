function [Rsquared_best, Rsquared_coeffs] = review_paper_plot(filename_dir,filename,Table,Align_data,Names,cutoff,cut,featureEND)
%   REVIEW_PAPER_PLOT runs the Linear Correlation experiment
%   -cut is the identity score of the lower bound for the segment in
%   quesiton
%   -cutoff is the row in the sorted Table that corresponds with cut
    
start=cutoff;
Features=[1:featureEND];
Predict_all=Table;
Align_all=Align_data(:,1);
Predict=Predict_all([start:end],:);
Align=Align_all(start:end);

numF=size(Features,2);
Rsquared_coeffs=zeros(numF,2);

for i=1:numF
    mdl=fitlm(Predict(:,Features(i)),Align);
    Rsquared_coeffs(i,2)=mdl.Rsquared.Ordinary;
    Rsquared_coeffs(i,1)=Features(i);
    i;
end


[~, I]=sort(Rsquared_coeffs(:,2));
Rsquared_coeffs=Rsquared_coeffs(I,:);
Features=[Rsquared_coeffs(numF,1) Rsquared_coeffs(numF-1,1) Rsquared_coeffs(numF-2,1) Rsquared_coeffs(numF-3,1) Rsquared_coeffs(numF-4,1)];

for k=1:5
    F=figure();
    plot(Align,Predict(:,Features(k)),'.');
    lsline;
    title(Names{Features(k)});
    xlabel('Identity Score');
    ylabel('Statistic');
    if cutoff~=1
        xlim([cut 1]);
    else
        xlim([.5 1]);
    end
    ylim([0 1]);
    set(gca,'fontsize',30);
    saveas(F, [filename_dir '/Fig/' num2str(cut*100) '-' num2str(k)] , 'epsc');
    
end

Rsquared_best=Rsquared_coeffs([numF numF-1 numF-2 numF-3 numF-4],1);
[~,J]=sort(Rsquared_coeffs(:,2),'descend');
Rsquared_coeffs=Rsquared_coeffs(J,:);
writefile_names(filename,Rsquared_coeffs,1,Names,'Linear');

disp('Top 5 statistcs:');
for i=1:5
    
fprintf('%s\n',Names{Features(i)});
end
fprintf('\n');



end