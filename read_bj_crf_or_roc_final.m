function []=read_bj_crf_or_roc_final(datamat,chNum,psychoname,testContrast,sampleContrast,animal,area,startEndTime,analysisType,excludeSessHighSSE,excludeOutliers,rootFolder,plotLeastSquares,modifyMinMax,useISI,calcPartial)
%Modified from read_blanco_V1_crf
%Written by Xing 06/03/13
%
%Reads from file containing average spike activity levels across sessions,
%for an individual channel at a time, as well as the two time points
%flanking each period. Performs contrast response function curve fitting
%with Naka Rushton function for each period for each session, generating as many
%curves for each session as there are time windows to analyse. 
%Also reads in data from file containing performance values, which is common
%to all channels, and has one set of 14 values per session. Performs
%psychometric curve fitting, generates 1 curve for each session.
%Calculates value of slopes of psychometric and CRF for
%each session, generates a graph combining slopes across sessions, for
%CRF data, and another graph for psychometric data. Plots a third
%graph of CRF slopes against psychometric slopes.
%Checks for correlations between
%1. CRF slope and time 2. Psychometric slope and time
%3. CRF against psychometric slope.
%Writes correlation coefficients to file, crf_correlation_coefs, as well as
%values of slopes and session numbers. Coefs: 1,1 CRF with time 1,2 Psycho with
%time, 1,3 CRF with Psycho. p-vals: 2,1 CRF with time 2,2 Psycho with
%time, 2,3 CRF with Psycho.

writeCoefs=1;
plotFig=1;
plotPsychoFig=0;
if strcmp(analysisType,'ROC_diff')||strcmp(analysisType,'ROC_diff2')
    plotPsychoFig=0;
end
excludeSessions=[26 50 306 312 316 322:328 342];
excludeSessions=[26 50 312 342 398 451];
% loadText=['load ',folder,'\allChROC',appendText,'.mat allChROC'];
% eval(loadText)
SSEcutoff=0.09;
slSigmaMultiple=3;
c50SigmaMultiple=3;
threshSigmaMultiple=3;
calculateTangent=1;
plotDiffC50_30=1;
            
appendText=['_',num2str(sampleContrast)];

if strcmp(analysisType,'ROC')||strcmp(analysisType,'NVP')||strcmp(analysisType,'ROC_diff2')||strcmp(analysisType,'ROC_zero_one')||strcmp(analysisType,'NVP_zero_one')
    manualCutoffMatText=['load F:\PL\ROC_mat_files\',animal,'\manual_cutoff.mat manualCutoff'];
    eval(manualCutoffMatText);
    if ~isempty(manualCutoff)
        ind=find(chNum==manualCutoff(:,1));
        if ~isempty(ind)
            manual_cutoff=manualCutoff(ind,2);
        else
            manual_cutoff=100;
        end
    else
        manual_cutoff=100;
    end
end

sessionSorted1=cell2mat(datamat(:,1))';
numsessions=length(datamat);
chSSE=zeros(length(sessionSorted1),2);
SSEMatFileName=[num2str(chNum),appendText,startEndTime,'_SSE_',area];
SSEMatFolder=fullfile('F:','PL',analysisType,animal,'SSE_mat_files');
if ~exist(SSEMatFolder,'dir')
    mkdir(SSEMatFolder);
end
SSEMatPath=fullfile(SSEMatFolder,SSEMatFileName);

if strcmp(analysisType,'ROC')||strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC_diff2')||strcmp(analysisType,'ROC_diff')||strcmp(analysisType,'ROC_zero_one')
    slC50Matname=[num2str(chNum),appendText,startEndTime,'_slC50_',area];
    slC50MatFolder=fullfile('F:','PL',analysisType,animal,'slope_C50_mat');
elseif strcmp(analysisType,'NVP')||strcmp(analysisType,'NVP_zero_one')
    slC50Matname=[num2str(chNum),appendText,startEndTime,'_nvpThreshold_',area];
    slC50MatFolder=fullfile('F:','PL',analysisType,animal,'nvpThreshold_mat');
end
if ~exist(slC50MatFolder,'dir')
    mkdir(slC50MatFolder);
end
slC50MatPathname=fullfile(slC50MatFolder,slC50Matname);

if plotFig==1&&~modifyMinMax
    [slopeNeuro,c50,diffc50,minRate,maxRate,chSSE,yLimData,threshold82lower,threshold82higher]=plot_CRF_or_ROC_across_sessions_final(animal,area,analysisType,datamat,chNum,numsessions,sessionSorted1,sampleContrast,testContrast,calculateTangent,plotDiffC50_30,excludeSessHighSSE,excludeOutliers,SSEMatPath,startEndTime,slC50MatPathname,slSigmaMultiple,c50SigmaMultiple,threshSigmaMultiple,rootFolder,plotLeastSquares,useISI);
end

% allChROC=[allChROC;appendROC];
% if excludeSessHighSSE==0
%     saveText=['save ',folder,'\allChROC.mat allChROC'];
% elseif excludeSessHighSSE==1
%     saveText=['save ',folder,'\allChROC_goodSSE.mat allChROC'];
% end
% eval(saveText)

if plotPsychoFig==1&&~strcmp(analysisType,'ROC_diff')
    plot_psycho_across_sessions(psychoname,sampleContrast,testContrast,excludeSessions,calculateTangent)
end

if modifyMinMax==1
    if isempty(chNum)
        chText='mean_across_channels';
    else
        chText=num2str(chNum);
    end
    %load previous minRate & maxRate values and calculate correct ones:
    if excludeSessHighSSE==0
        coefMatname=[chText,appendText,startEndTime,'_',analysisType,'_coefs_',area];
    elseif excludeSessHighSSE==1
        if excludeOutliers==0
            coefMatname=[chText,appendText,startEndTime,'_',analysisType,'_coefs_',area,'_goodSSE'];
        elseif excludeOutliers==1
            coefMatname=[chText,appendText,startEndTime,'_',analysisType,'_coefs_',area,'_goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
        end
    end
    subFolderName=[analysisType,'_coef_mat'];
    coefMatFolder=fullfile('F:','PL',analysisType,animal,subFolderName);
    coefMatPathname=fullfile(coefMatFolder,coefMatname);
    if ~exist(coefMatFolder,'dir')
        mkdir(coefMatFolder)
    end
    if plotDiffC50_30==1
        loadText=['load ',coefMatPathname,'.mat coefficients sessionSorted1 slopeNeuro sessionSorted2 slopePsycho matchPsycho c50 diffc50 minRate maxRate sessionSorted1 threshold82higher'];
    else
        loadText=['load ',coefMatPathname,'.mat coefficients sessionSorted1 slopeNeuro sessionSorted2 slopePsycho matchPsycho c50 minRate maxRate sessionSorted1 threshold82higher'];
    end
    eval(loadText)
    incorrectMaxRate=maxRate;
    incorrectMinRate=minRate;
    maxRate=1-incorrectMinRate;
    minRate=1-incorrectMinRate-incorrectMaxRate;
end
if strcmp(analysisType,'ROC')||strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC_zero_one')
    example_ch_54=0;
    plot_neurometric_coefs_final(animal,area,chNum,appendText,startEndTime,slopeNeuro,c50,plotDiffC50_30,diffc50,minRate,maxRate,sessionSorted1,analysisType,example_ch_54,excludeSessHighSSE,excludeOutliers,writeCoefs,slSigmaMultiple,c50SigmaMultiple,calcPartial,threshold82higher,threshold82lower)
elseif strcmp(analysisType,'NVP_zero_one')
    plot_nvp_threshold_coefs_final(animal,area,chNum,appendText,startEndTime,threshold82lower,threshold82higher,sessionSorted1,analysisType,excludeSessHighSSE,excludeOutliers,writeCoefs,threshSigmaMultiple,useISI)
elseif strcmp(analysisType,'NVP')
    plot_nvp_threshold_coefs_final(animal,area,chNum,appendText,startEndTime,[],threshold82higher,sessionSorted1,analysisType,excludeSessHighSSE,excludeOutliers,writeCoefs,threshSigmaMultiple,useISI)
end
% pause
if ~strcmp(analysisType,'ROC_diff2')
    close all
end