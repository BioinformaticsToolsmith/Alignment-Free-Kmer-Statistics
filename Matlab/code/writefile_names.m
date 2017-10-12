function CellTable = writefile_names(filename,Data,col,Names,KNNorLinear)
%converts a Data set to a file
%returns a cell array of the table
%col is the column where the statistic numbers are converted to names

rows=size(Data,1);
CellTable=num2cell(Data);

for i=1:rows
    CellTable{i,1}=Names{Data(i,col)};
end

if strcmp(KNNorLinear,'KNN')
    T=cell2table(CellTable,'VariableNames',{'Statistics' 'Correct_Neighbors'});
    writetable(T,filename,'Delimiter','\t');
elseif strcmp(KNNorLinear,'Linear')
    T=cell2table(CellTable,'VariableNames',{'Statistics' 'R_squared_coefficient'});
    writetable(T,filename,'Delimiter','\t');
else
    T=cell2table(CellTable,'VariableNames',{'Statistics' 'Time_seconds'});
    writetable(T,filename,'Delimiter','\t');
end

end