function EVALUATE(path)

pack;
disp(['Turning off warnings']);
warning ('off','all');

%   EVALUATE is the overall driver program
%
%
Hist1_path=[path '/Hist/hist1.out'];         %these histograms are required
Hist2_path=[path '/Hist/hist2.out'];         % for certains statistics
Hist1=dlmread(Hist1_path,'!',0,1);
Hist2=dlmread(Hist2_path,'!',0,1);
sizes=sum(Hist1,2);
ave_size=mean(sizes);
k=round(log(ave_size)/log(4));            %selects an appropriate k-mer size
disp(k);
switch(k)
    case 1
        Histk_path=[path '/Hist/hist1.out'];
        Histk=dlmread(Histk_path,'!',0,1);
        Histk_=[];
        Histk__=[];
        Histkp_path=[path '/Hist/hist2.out'];
        Histkp=dlmread(Histkp_path,'!',0,1);
        Histkpp_path=[path '/Hist/hist3.out'];
        Histkpp=dlmread(Histkpp_path,'!',0,1);

    case 2
        Histk_path=[path '/Hist/hist2.out'];
        Histk=dlmread(Histk_path,'!',0,1);
        Histk_path2=[path '/Hist/hist1.out'];
        Histk_=dlmread(Histk_path2,'!',0,1);
        Histk__=[];
        Histkp_path=[path '/Hist/hist3.out'];
        Histkp=dlmread(Histkp_path,'!',0,1);
        Histkpp_path=[path '/Hist/hist4.out'];
        Histkpp=dlmread(Histkpp_path,'!',0,1);

    case 3
        Histk_path=[path '/Hist/hist3.out'];
        Histk=dlmread(Histk_path,'!',0,1);
        Histk_path2=[path '/Hist/hist2.out'];
        Histk_=dlmread(Histk_path2,'!',0,1);
        Histk_path3=[path '/Hist/hist1.out'];
        Histk__=dlmread(Histk_path3,'!',0,1);
        Histkp_path=[path '/Hist/hist4.out'];
        Histkp=dlmread(Histkp_path,'!',0,1);
        Histkpp_path=[path '/Hist/hist5.out'];
        Histkpp=dlmread(Histkpp_path,'!',0,1);

    case 4
        Histk_path=[path '/Hist/hist4.out'];
        Histk=dlmread(Histk_path,'!',0,1);
        Histk_path2=[path '/Hist/hist3.out'];
        Histk_=dlmread(Histk_path2,'!',0,1);
        Histk_path3=[path '/Hist/hist2.out'];
        Histk__=dlmread(Histk_path3,'!',0,1);
        Histkp_path=[path '/Hist/hist5.out'];
        Histkp=dlmread(Histkp_path,'!',0,1);
        Histkpp_path=[path '/Hist/hist6.out'];
        Histkpp=dlmread(Histkpp_path,'!',0,1);

    case 5
        Histk_path=[path '/Hist/hist5.out'];
        Histk=dlmread(Histk_path,'!',0,1);
        Histk_path2=[path '/Hist/hist4.out'];
        Histk_=dlmread(Histk_path2,'!',0,1);
        Histk_path3=[path '/Hist/hist3.out'];
        Histk__=dlmread(Histk_path3,'!',0,1);
        Histkp_path=[path '/Hist/hist6.out'];
        Histkp=dlmread(Histkp_path,'!',0,1);
        Histkpp_path=[path '/Hist/hist7.out'];
        Histkpp=dlmread(Histkpp_path,'!',0,1);

    case 6
        Histk_path=[path '/Hist/hist6.out'];
        Histk=dlmread(Histk_path,'!',0,1);
        Histk_path2=[path '/Hist/hist5.out'];
        Histk_=dlmread(Histk_path2,'!',0,1);
        Histk_path3=[path '/Hist/hist4.out'];
        Histk__=dlmread(Histk_path3,'!',0,1);
        Histkp_path=[path '/Hist/hist7.out'];
        Histkp=dlmread(Histkp_path,'!',0,1);
        Histkpp_path=[path '/Hist/hist8.out'];
        Histkpp=dlmread(Histkpp_path,'!',0,1);
    case 7
        Histk_path=[path '/Hist/hist7.out'];
        Histk=dlmread(Histk_path,'!',0,1);
        Histk_path2=[path '/Hist/hist6.out'];
        Histk_=dlmread(Histk_path2,'!',0,1);
        Histk_path3=[path '/Hist/hist5.out'];
        Histk__=dlmread(Histk_path3,'!',0,1);
        Histkp_path=[path '/Hist/hist8.out'];
        Histkp=dlmread(Histkp_path,'!',0,1);
        Histkpp_path=[path '/Hist/hist9.out'];
        Histkpp=dlmread(Histkpp_path,'!',0,1);

    case 8
        Histk_path=[path '/Hist/hist8.out'];
        Histk=dlmread(Histk_path,'!',0,1);
        Histk_path2=[path '/Hist/hist7.out'];
        Histk_=dlmread(Histk_path2,'!',0,1);
        Histk_path3=[path '/Hist/hist6.out'];
        Histk__=dlmread(Histk_path3,'!',0,1);
        Histkp_path=[path '/Hist/hist9.out'];
        Histkp=dlmread(Histkp_path,'!',0,1);
        Histkpp_path=[path '/Hist/hist10.out'];
        Histkpp=dlmread(Histkpp_path,'!',0,1);

    case 9
        Histk_path=[path '/Hist/hist9.out'];
        Histk=dlmread(Histk_path,'!',0,1);
        Histk_path2=[path '/Hist/hist8.out'];
        Histk_=dlmread(Histk_path2,'!',0,1);
        Histk_path3=[path '/Hist/hist7.out'];
        Histk__=dlmread(Histk_path3,'!',0,1);
        Histkp_path=[path '/Hist/hist10.out'];
        Histkp=dlmread(Histkp_path,'!',0,1);
        Histkpp_path=[path '/Hist/hist11.out'];
        Histkpp=dlmread(Histkpp_path,'!',0,1);

    case 10
        Histk_path=[path '/Hist/hist10.out'];
        Histk=dlmread(Histk_path,'!',0,1);
        Histk_path2=[path '/Hist/hist9.out'];
        Histk_=dlmread(Histk_path2,'!',0,1);
        Histk_path3=[path '/Hist/hist8.out'];
        Histk__=dlmread(Histk_path3,'!',0,1);
        Histkp_path=[path '/Hist/hist11.out'];
        Histkp=dlmread(Histkp_path,'!',0,1);
        Histkpp_path=[path '/Hist/hist12.out'];
        Histkpp=dlmread(Histkpp_path,'!',0,1);

    otherwise
        error('bad k value');
end

Align_path=[path '/Align/align.out'];
Align=dlmread(Align_path,'!',0,3);

numSeqs=size(Hist1,1);
Align_data=addSeqNums(Align,1,numSeqs);  %add seq numbers to the alignment identity score data

% disp(['This is the path for the files: ' path]);
% disp('running Matlab evaluation');

REVIEW_PAPER(Align_data,Histk,Histk_,Histk__,Histkp,Histkpp,Hist1,Hist2,path);

disp('Done');
end