function Names = make_feature_names(List)
%   MAKE_FEATURE_NAMES returns a cell of all the given paired-statistic
%   names
%   -List is a cell of each of the single statistic names


N=size(List,2);
Names=cell(N*2,1);
for n=1:N
   Names{n}=List{n} ;
end
a=1;
for i=N+1:2*N
   Names{i}=strcat(List{a},'^2');
   a=a+1;
end
stop=2*N;
for j=1:stop
    f1=Names{j};
    for m=j+1:stop
        f2=Names{m};
        result1=strcat(f1,' * ');
        result=strcat(result1, f2);
        Names=[Names; result];
    end
end


end