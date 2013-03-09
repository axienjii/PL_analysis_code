function []=read_bj_V4_1_roc(rocname,psychoname,testContrast,sampleContrast,skipSessions,chNum,area,folder,appendText,animal,sessions)
%created on 05/02/13 for GitHub version control.
%Modified from read_jack_V4_1_roc.
%Combines ROC values across channels into mat file, allChROC.mat.
%Reads from file containing ROC values (14 per session) across sessions,
%for an individual channel at a time. Performs neurometric curve fitting 
%with Weibull function for each session, generating 1 curve for each session. 
%If SSE for fitted curve is over a certain value, exclude that session.
%Also reads in data from file containing performance values, which is common 
%to all channels, and has one set of 14 values per session. Performs 
%psychometric curve fitting, generates 1 curve for each session. 
%Calculates value of slopes of psychometric and neurometric functions for 
%each session, generates a graph combining slopes across sessions, for 
%neurometric data, and another graph for psychometric data. Plots a third 
%graph of neurometric slopes against psychometric slopes. 
%Checks for correlations between 
%1. Neurometric slope and time 2. C50 slope and time 
%3. Neurometric against psychometric slope. 
%Writes correlation coefficients to file, roc_correlation_coefs, as well as 
%values of slopes and session numbers. Coefs: 1,1 Neuro with time 1,2 C50 with
%time, 1,3 Neuro with Psycho. p-vals: 2,1 Neuro with time 2,2 C50 with
%time, 2,3 Neuro with Psycho.

writeCoefs=1;
plotFig=1;
plotPsychoFig=1;
excludeSessions=[26 50 306 312 316 322:328 342];
% loadText=['load ',folder,'\allChROC',appendText,'.mat allChROC'];
% eval(loadText)
rocMatFileName=[rocname,'_roc.mat'];
excludeSessHighSSE=1;%set to 1 to exclude sessions with poor Weibull fit
SSEcutoff=0.09;
excludeOutliers=1;
slSigmaMultiple=3;
c50SigmaMultiple=3;
calculateTangent=1;
plotDiffC50_30=1;
            
manualCutoffMatText=['load F:\PL\ROC_mat_files\',animal,'\manual_cutoff.mat manualCutoff'];
eval(manualCutoffMatText);
ind=find(chNum==manualCutoff(:,1));
if ~isempty(ind)
    manual_cutoff=manualCutoff(ind,2);
else
    manual_cutoff=100;
end

loadMatText=['load ',rocMatFileName];
eval(loadMatText)
goodCount=1;
for rowNum=1:size(rocValsMat,1)
    ind=find(skipSessions(:,1)==chNum);
    skip=0;
    if ~isempty(ind)
        if find(skipSessions(ind,2:end)==rocValsMat{rowNum,2})
            skip=1;
        end
    end
    if skip==0
        sessionRow=rocValsMat{rowNum,2};
        if ~isempty(find(sessionRow==sessions))
            VALUES(goodCount,:)=rocValsMat{rowNum,3};
            session1(goodCount)=rocValsMat{rowNum,2};
            goodCount=goodCount+1;
        end
    end
end
    
for excludeCount=1:length(excludeSessions)
    ind=(session1==excludeSessions(excludeCount));
    VALUES=VALUES(~ind,:);
    session1=session1(~ind);
end

    [sessionSorted1 index]=sort(session1);
%     sessionSorted1=session1;%sessions are in order so no need to sort
    numsessions=size(VALUES,1);
    for i=1:numsessions
%         VALUESsorted(i,1:16)=VALUES(index(i),1:16);
        VALUESsorted(i,:)=VALUES(index(i),:);
    end
    VALUES=VALUESsorted;

SSEMatFileName=[num2str(chNum),'_',num2str(sampleContrast),'_SSE'];
SSEMatFolder=fullfile('F:','PL','SSE_mat_files',animal);
SSEMatPath=fullfile(SSEMatFolder,SSEMatFileName);
slC50Matname=[num2str(chNum),appendText,'_slC50'];
slC50MatPathname=fullfile('F:','PL','slope_C50_mat',animal,slC50Matname);
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
    [slopeNeuro,c50,diffc50,minRate,maxRate]=plot_CRF_or_ROC_across_sessions(animal,'ROC',ROCmat,chNum,numsessions,sessionSorted1,sampleContrast,testContrast,calculateTangent,plotDiffC50_30,excludeSessHighSSE,excludeOutliers,SSEMatPath,startEndTime);
end

if plotPsychoFig==1
    plot_psycho_across_sessions(psychoname,sampleContrast,testContrast,excludeSessions,calculateTangent)
end

if plotDiffC50_30==1
    fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.5, 0.9]);
    set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21 28.41]);
else
    fighandle2=figure('Color',[1,1,1],'Units','Normalized','Position',[0.14, 0.46, 0.8, 0.4]);
    set(fighandle2, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
    orient landscape
end

example_ch_54=0;
analysisTypeText='ROC';
plot_neurometric_coefs(animal,area,chNum,appendText,startEndTime,slopeNeuro,c50,plotDiffC50_30,diffc50,minRate,maxRate,sessionSorted1,analysisTypeText,example_ch_54,excludeSessHighSSE,excludeOutliers,writeCoefs)

% pause
close all
