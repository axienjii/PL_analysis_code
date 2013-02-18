function [SE_TimeStamps,missing]=open_nse_file(nsePath)
%Written 31/01/13 by Xing
missing=0;
SE_TimeStamps=[];ScNumbers=[];CellNumbers=[];Params=[];DataPoints=[];NlxHeader=[];
if exist(nsePath,'file')
    [SE_TimeStamps,ScNumbers,CellNumbers,Params,DataPoints,NlxHeader]=Nlx2MatSpike(nsePath,[1 1 1 1 1],1,1,1);
else
    missing=1;
    load('F:\PL\missingNSEfile.mat')
    missingNSEfile=[missingNSEfile;{nsePath}];
    save F:\PL\missingNSEfile.mat missingNSEfile
end