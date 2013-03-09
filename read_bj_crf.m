function []=read_bj_crf(CRFmat,chNum,psychoname,testContrast,sampleContrast,animal,area,startEndTime)
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
plotPsychoFig=1;
excludeSessions=[26 50 306 312 316 322:328 342];
% loadText=['load ',folder,'\allChROC',appendText,'.mat allChROC'];
% eval(loadText)
excludeSessHighSSE=0;%set to 1 to exclude sessions with poor Weibull fit
SSEcutoff=0.09;
excludeOutliers=0;
slSigmaMultiple=3;
c50SigmaMultiple=3;
calculateTangent=1;
plotDiffC50_30=1;
            
appendText=['_',num2str(sampleContrast)];

% manualCutoffMatText=['load F:\PL\ROC_mat_files\',animal,'\manual_cutoff.mat manualCutoff'];
% eval(manualCutoffMatText);
% ind=find(chNum==manualCutoff(:,1));
% if ~isempty(ind)
%     manual_cutoff=manualCutoff(ind,2);
% else
%     manual_cutoff=100;
% end

sessionSorted1=cell2mat(CRFmat(:,1));
numsessions=length(CRFmat);
chSSE=zeros(length(sessionSorted1),2);
fighandle1=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);

SSEMatFileName=[num2str(chNum),'_',appendText,startEndTime,'_SSE'];
SSEMatFolder=fullfile('F:','PL','CRF',animal,'SSE_mat_files');
if ~exist(SSEMatFolder,'dir')
    mkdir(SSEMatFolder);
end
SSEMatPath=fullfile(SSEMatFolder,SSEMatFileName);
slC50Matname=[num2str(chNum),'_',appendText,startEndTime,'_slC50'];
slC50MatFolder=fullfile('F:','PL','CRF',animal,'slope_C50_mat');
if ~exist(slC50MatFolder,'dir')
    mkdir(slC50MatFolder);
end
slC50MatPathname=fullfile(slC50MatFolder,slC50Matname);
if excludeSessHighSSE==1
    loadText=['load ',SSEMatPath,' chSSE'];
    eval(loadText);
    SSEcutoff=mean(chSSE(:,2))+std(chSSE(:,2));
    ind=(chSSE(:,2)<SSEcutoff);
    sessionSorted1=sessionSorted1(ind);
    VALUES=VALUES(ind,:);
    numsessions=length(sessionSorted1);
    if excludeOutliers==1   %within reduced pool of sessions with good SSE, further examine slope and C50 values for outliers
        loadText=['load ',slC50MatPathname,'.mat sessionSorted1 slopeNeuro c50'];
        eval(loadText)    
        c50=real(c50);   
        slsigma=std(slopeNeuro);
        c50sigma=std(c50);
        sloutliers=abs((slopeNeuro-mean(slopeNeuro)))>slSigmaMultiple*slsigma;
        c50outliers=abs((c50-mean(c50)))>c50SigmaMultiple*c50sigma;
        c50outliersHighcut=c50>manual_cutoff;%manual exclusion for obvious outliers
        c50outliersLowcut=c50<0;%manual exclusion for obvious outliers
%         if sum(c50outliersHighcut)>0
%             pause
%         end
        ind=sloutliers+c50outliers+c50outliersHighcut+c50outliersLowcut;%find sessions where slope and/or C50 values are outliers (union)
        sessionSorted1=sessionSorted1(~ind);%keep sessions that do not have outliers
        VALUES=VALUES(~ind,:);
        numsessions=length(sessionSorted1);
        slopeNeuro=slopeNeuro(~ind);       
        c50=c50(~ind);       
    end
end

if plotFig==1
    [slopeNeuro,c50,diffc50,minRate,maxRate]=plot_CRF_or_ROC_across_sessions(animal,area,'CRF',CRFmat,chNum,numsessions,sessionSorted1,sampleContrast,testContrast,calculateTangent,plotDiffC50_30,excludeSessHighSSE,excludeOutliers,SSEMatPath,startEndTime);
end

% allChROC=[allChROC;appendROC];
% if excludeSessHighSSE==0
%     saveText=['save ',folder,'\allChROC.mat allChROC'];
% elseif excludeSessHighSSE==1
%     saveText=['save ',folder,'\allChROC_goodSSE.mat allChROC'];
% end
% eval(saveText)

if plotPsychoFig==1
    plot_psycho_across_sessions(psychoname,sampleContrast,testContrast,excludeSessions,calculateTangent)
end

example_ch_54=0;
analysisTypeText='CRF';
plot_neurometric_coefs(animal,area,chNum,appendText,startEndTime,slopeNeuro,c50,plotDiffC50_30,diffc50,minRate,maxRate,sessionSorted1,analysisTypeText,example_ch_54,excludeSessHighSSE,excludeOutliers,writeCoefs)

% pause
close all