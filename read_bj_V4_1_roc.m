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
            
xvals=testContrast(1):1:testContrast(end);

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
X0=[2 30 0.2 0.1];
if chNum==13||chNum==24||chNum==42||chNum==52%stimulus-induced suppression, value of slope is negative
    X0=[-2 30 0.2 0.1];
end
chSSE=zeros(length(sessionSorted1),2);
options = optimset('Display','off','MaxFunEvals',10^8,'MaxIter',10^6,'TolFun',1.0E-6,'TolX',1.0E-6);
fighandle1=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);

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

for i=1:numsessions
    rocvals=VALUES(i,1:length(testContrast));%for channel of interest, for each session
    subplot(ceil(numsessions/5),5,i);
%     subplot(3,numsessions,i);%V1 data
    if sum(rocvals(1:3))<sum(rocvals(end-2:end))||chNum==24&&sessionSorted1(1,i)==335||chNum==24&&sessionSorted1(1,i)==336||chNum==13
        X0=[2 30 0.2 0.1];
    elseif sum(rocvals(1:3))>sum(rocvals(end-2:end))
        X0=[-2 30 0.2 0.1];
    end
    X=fminsearch(@fit_weibull,X0,[],testContrast,rocvals,[],'least_square',[1 1 0 0],[20 100 0 0],[1 1 0 0],[-20 0 0 0]);
    if calculateTangent==0
        slopeNeuro(1,i)=X(1);
    elseif calculateTangent==1
        slopeNeuro(1,i)=X(1)*X(3)*exp(-(30/X(2))^X(1) )*30^(X(1)-1)*(1/X(2))^X(1);
    end
    c50(1,i)=real(X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1)));
%     [X,fval]=fminsearch('weib_sim_min_max',X0,options,testContrast,rocvals);
%     slopeNeuro(1,i)=X(2);      
%     c50(1,i)=X(1).*(-log(0.5)).^(1/X(2));
%     yvals=max(rocvals)-(max(rocvals)-min(rocvals)).*exp(-((xvals./X(1)).^X(2)));
    plot(testContrast,rocvals,'ok');
    hold on
    if plotDiffC50_30==1
        diffc50(1,i)=abs(c50(1,i)-30);
    end
        line(c50(1,i),0:0.01:1,'Color','r');
    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
    fitted_yvals=1-X(4)-X(3).*exp(-(testContrast./X(2)).^X(1));
    residuals = rocvals-fitted_yvals;
	sseROC=sum(residuals.^2);
    chSSE(i,:)=[sessionSorted1(i) sseROC];

    plot(xvals,yvals,'r');
    line(sampleContrast,0:0.01:1);
    ylim([0,max(rocvals)]);
% ylim([0 1]);
	xlim([0 max(testContrast)]);
    subplottitle=num2str(sessionSorted1(1,i));
    title(subplottitle);
    if i==1
        ptext=sprintf('%s',rocname);
        orient landscape
        yLimVals=get(gca,'YLim');
        text('Position',[-10 yLimVals(2)+0.2],'FontSize',9,'String',ptext);
    end
    ptext=sprintf('SSE: %.3f',sseROC);
    text('Position',[testContrast(end)/2 yLimVals(1)+0.05],'FontSize',9,'String',ptext);
end
if excludeSessHighSSE==0
    saveText=['save ',SSEMatPath,' chSSE'];
    eval(saveText);
    rocImagename=[num2str(chNum),appendText];
elseif excludeSessHighSSE==1
    if excludeOutliers==0
        rocImagename=[num2str(chNum),appendText,'_goodSSE'];
        saveText=['save ',slC50MatPathname,'.mat sessionSorted1 slopeNeuro c50'];
        eval(saveText)
    elseif excludeOutliers==1
        rocImagename=[num2str(chNum),appendText,'_goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
    end
end
rocPathname=fullfile('F:','PL','neurometric_ROC',animal,rocImagename);
printtext=sprintf('print -dpng %s.png',rocPathname);
eval(printtext);

appendROC(1:numsessions,1)=chNum;
appendROC(1:numsessions,2)=sessionSorted1;
appendROC(1:numsessions,3:length(testContrast)+2)=VALUES;
if plotDiffC50_30==1
    appendROC=[appendROC slopeNeuro' c50' diffc50'];
else
    appendROC=[appendROC slopeNeuro' c50'];
end
% allChROC=[allChROC;appendROC];
% if excludeSessHighSSE==0
%     saveText=['save ',folder,'\allChROC.mat allChROC'];
% elseif excludeSessHighSSE==1
%     saveText=['save ',folder,'\allChROC_goodSSE.mat allChROC'];
% end
% eval(saveText)

fid=fopen(psychoname,'r');%read file and calculate number of sessions
[A,count1]=fscanf (fid,'%s ', inf);
numsessions=count1/(length(testContrast)+5)
fclose(fid);

count=1;
session2=zeros(1,numsessions);
slopePsycho=zeros(1,numsessions);

% X0=[30 2];
X0=[2 30 0.2 0.1];
fid=fopen(psychoname,'r');
while (count<=numsessions)
    [ID,N]=fscanf (fid,'%d', 1);%and session ID
    session2(1,count)=ID;
    [values,count2]=fscanf(fid,'%f', [1,length(testContrast)+4]);
    VALUES(count,1:length(testContrast)+4)=values;
    count=count+1;
end;

for excludeCount=1:length(excludeSessions)
    ind=(session2==excludeSessions(excludeCount));
    VALUES=VALUES(~ind,:);
    session2=session2(~ind);
end

[sessionSorted2 index]=sort(session2);
VALUESsorted=[];
numsessions=length(session2);
for i=1:numsessions
        VALUESsorted(i,:)=VALUES(index(i),:);
end
VALUES=VALUESsorted;

figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
subplotInd=0;
for i=1:numsessions
    subplotInd=subplotInd+1;
    perfvals=VALUES(i,1:length(testContrast));%psychometric performance for each session
    X=fminsearch(@fit_weibull,X0,[],testContrast,perfvals,[],'least_square',[0 0 0 0],[],[0 0 0 0],[]);
    %         [X,fval]=fminsearch('weib_sim_min_max',X0,options,testContrast,perfvals);
    %         slopePsycho(1,i)=X(2);
    %         PSE=X(1).*(-log(0.5)).^(1/X(2));
    %         yvals=max(perfvals)-(max(perfvals)-min(perfvals))*exp(-((xval
    %         s/X(1)).^X(2)));
    if calculateTangent==0
        slopePsycho(1,i)=X(1);
    elseif calculateTangent==1
        slopePsycho(1,i)=X(1)*X(3)*exp(-(30/X(2))^X(1) )*30^(X(1)-1)*(1/X(2))^X(1);
    end
    PSE=X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1));
    subplot(ceil((numsessions)/5),5,subplotInd);
    plot(testContrast,perfvals,'ok');
    hold on
    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
    plot(xvals,yvals,'r');
    line(PSE,0:0.01:1,'Color','r');
    line(sampleContrast,0:0.01:1);
    ylim([0,max(perfvals)]);
    xlim([0 max(testContrast)]);
    subplottitle=num2str(sessionSorted2(1,i));
    title(subplottitle);
    if i==1
        ptext=sprintf('%s',psychoname);
        orient landscape
        yLimVals=get(gca,'YLim');
        text('Position',[-10 yLimVals(2)+0.2],'FontSize',9,'String',ptext);
    end
end
printtext=sprintf('print -dpng %s',psychoname);
eval(printtext);

if plotDiffC50_30==1
    fig=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.04, 0.5, 0.9]);
    set(fig,'PaperUnits','centimeters','PaperType','A4','PaperOrientation', 'portrait', 'PaperPosition', [0.63452 0.63452 21 28.41]);
else
    fighandle2=figure('Color',[1,1,1],'Units','Normalized','Position',[0.14, 0.46, 0.8, 0.4]);
    set(fighandle2, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
    orient landscape
end
%plot slope vals against ordinal, not actual, session numbers:
% test=subplot(1,3,1);
% plot(1:length(sessionSorted1),slopeNeuro,'ok');
% test=subplot(1,3,2);
% plot(1:length(sessionSorted2),slopePsycho,'ok');
if plotDiffC50_30==1
    subplot(2,2,1);
else
    subplot(1,3,1);
end
example_ch_54=0;
if example_ch_54==1
    xvals=0:20;
    plot(xvals,slopeNeuro,'sb','MarkerFaceColor','b');hold on
    brob = robustfit(xvals,slopeNeuro);
    plot(xvals,brob(1)+brob(2)*xvals,'k','LineWidth',4);hold on
    title('neurometric slope vs time')
    if plotDiffC50_30==1
        subplot(1,3,3);
        plot(xvals,c50,'sb');
        title('c50 vs time')
        subplot(2,2,4);
        plot(xvals,diffc50,'sb');
        title('abs(c50-30) vs time')
    else
        subplot(1,3,2);
        plot(xvals,c50,'sb');
        title('c50 vs time')
    end
    if plotDiffC50_30==1
        subplot(2,2,2);
    else
        subplot(1,3,3);
    end
    matchPsycho=zeros(1,length(slopeNeuro));%find slopePsycho for sessions where data from slopeNeuro is available
    for i=1:length(slopeNeuro)
        ind=find(sessionSorted1(i)==sessionSorted2);
        matchPsycho(i)=slopePsycho(ind(1));%if 2 cells present for 1 channel, simply duplicates value of slopePsycho for that session
    end
    plot(matchPsycho,slopeNeuro,'sb');
    title('neurometric slope vs psychometric slope')
else
%     plot(sessionSorted1,slopeNeuro,'ok');
    plot(1:length(sessionSorted1),slopeNeuro,'ok');
    title('neurometric slope vs time')
    if plotDiffC50_30==1
        subplot(2,2,3);
        plot(1:length(c50),c50,'ok');
        title('c50 vs time')
        subplot(2,2,4);
        plot(1:length(sessionSorted1),diffc50,'ok');
        title('abs(c50-30) vs time')
    else
        subplot(1,3,2);
        plot(1:length(sessionSorted1),c50,'ok');
        title('c50 vs time')
    end
%     plot(sessionSorted1,c50,'ok');    
    if plotDiffC50_30==1
        subplot(2,2,2);
    else
        subplot(1,3,3);
    end
    matchPsycho=zeros(1,length(slopeNeuro));%find slopePsycho for sessions where data from slopeNeuro is available
    for i=1:length(slopeNeuro)
        ind=find(sessionSorted1(i)==sessionSorted2);
        matchPsycho(i)=slopePsycho(ind(1));%if 2 cells present for 1 channel, simply duplicates value of slopePsycho for that session
    end
    plot(matchPsycho,slopeNeuro,'ok');
    title('neurometric slope vs psychometric slope')
end
% slopeNeuro
% slopePsycho
% matchPsycho
% sessionSorted1
% sessionSorted2
% session1
% session2
a=[sessionSorted1' slopeNeuro'];
c50=real(c50);
b=[sessionSorted1' c50'];
c=[slopeNeuro' matchPsycho'];
[coefficients1 p1]=corrcoef(a); 
[coefficients2 p2]=corrcoef(b); 
[coefficients3 p3]=corrcoef(c);
if plotDiffC50_30==1
    diffc50=real(diffc50);
    d=[sessionSorted1' diffc50'];
    [coefficients4 p4]=corrcoef(d);
end
if length(coefficients1)>1
if plotDiffC50_30==1
    coefficients(1,1:4)=[coefficients1(2) coefficients2(2) coefficients3(2) coefficients4(2)];
    coefficients(2,1:4)=[p1(2) p2(2) p3(2) p4(2)]
    subplot(2,2,1);
    ptext=sprintf('r= %f  p= %f',coefficients(1,1),coefficients(2,1));
    yLimVals=get(gca,'YLim');
    xLimVals=get(gca,'XLim');
    text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
    subplot(2,2,3);
    ptext=sprintf('r= %f  p= %f',coefficients(1,2),coefficients(2,2));
    yLimVals=get(gca,'YLim');
    xLimVals=get(gca,'XLim');
    text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
    subplot(2,2,2);
    ptext=sprintf('r= %f  p= %f',coefficients(1,3),coefficients(2,3));
    yLimVals=get(gca,'YLim');
    xLimVals=get(gca,'XLim');
    text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);  
    subplot(2,2,4);
    ptext=sprintf('r= %f  p= %f',coefficients(1,4),coefficients(2,4));
    yLimVals=get(gca,'YLim');
    xLimVals=get(gca,'XLim');
    text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);   
else
    coefficients(1,1:3)=[coefficients1(2) coefficients2(2) coefficients3(2)];
    coefficients(2,1:3)=[p1(2) p2(2) p3(2)]
    subplot(1,3,1);
    ptext=sprintf('r= %f  p= %f',coefficients(1,1),coefficients(2,1));
    yLimVals=get(gca,'YLim');
    xLimVals=get(gca,'XLim');
    text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
    subplot(1,3,2);
    ptext=sprintf('r= %f  p= %f',coefficients(1,2),coefficients(2,2));
    yLimVals=get(gca,'YLim');
    xLimVals=get(gca,'XLim');
    text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
    subplot(1,3,3);
    ptext=sprintf('r= %f  p= %f',coefficients(1,3),coefficients(2,3));
    yLimVals=get(gca,'YLim');
    xLimVals=get(gca,'XLim');
    text('Position',[xLimVals(1) yLimVals(1)-0.1*(yLimVals(2)-yLimVals(1))],'FontSize',9,'String',ptext);
end
    set(gcf,'PaperPositionMode','auto')
    if excludeSessHighSSE==0
        rocCoefImagename=[num2str(chNum),appendText,'_roc_coefs'];
    elseif excludeSessHighSSE==1
        if excludeOutliers==0
            rocCoefImagename=[num2str(chNum),appendText,'_roc_coefs_goodSSE'];
        elseif excludeOutliers==1
            rocCoefImagename=[num2str(chNum),appendText,'_roc_coefs_goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
        end
    end
    rocCoefPathname=fullfile('F:','PL','ROC_coef_images',animal,rocCoefImagename);
    printtext=sprintf('print -dpng %s.png',rocCoefPathname);
    eval(printtext);
    if writeCoefs==1
        %write correlation coefficients, session numbers, and slope values
        %to file:
        if excludeSessHighSSE==0
            rocCoefMatname=[num2str(chNum),appendText,'_roc_coefs_',area];
        elseif excludeSessHighSSE==1
            if excludeOutliers==0
                rocCoefMatname=[num2str(chNum),appendText,'_roc_coefs_',area,'_goodSSE'];
            elseif excludeOutliers==1
                rocCoefMatname=[num2str(chNum),appendText,'_roc_coefs_',area,'_goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
            end
        end
        rocCoefMatPathname=fullfile('F:','PL','ROC_coef_mat',animal,rocCoefMatname);
        saveText=['save ',rocCoefMatPathname,'.mat coefficients sessionSorted1 slopeNeuro sessionSorted2 slopePsycho matchPsycho c50 diffc50'];
        eval(saveText)
    end
end

% pause
close all
