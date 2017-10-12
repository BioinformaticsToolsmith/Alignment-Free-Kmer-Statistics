function T = addSeqNums(Table, start, num)
%   ADDSEQNUMS adds the sequence numbers to a given alignment table
%   -Table: the alignment identity scores in a column
%   -start: the row number in the table to start adding seq numbers, usually 1 
%   -num: the number of sequences
T_=zeros(size(Table,2),2);
count=1;
for n=start:(start+num-1)
   for m=n:(start+num-1)
      T_(count,1)=n;
      T_(count,2)=m;
      count=count+1;
   end
    
end
T=[Table T_];

end