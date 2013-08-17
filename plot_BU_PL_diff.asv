function plot_BU_PL_diff(PL,BU)
%Written by Xing on 13/11/12. Plots ROC curves for PL stimuli, during
%performance of contrast discrimination task (attend to PL stmiuli in RF,
%blue), and during performance of orientation discrimination task (attend
%away from RF, red).
%Calculates difference between ROC values of the two attention tasks, and
%difference index, for each channel and session.

calculateTangent=1;
readData=0;
samePlot=1;
plotIndivChs=0;
areas=['v1_3';'v4_3'];
areaTexts=[{'V1'} {'V4'}];
for areaCount=1:2
    area=areas(areaCount,:);
    folderPrint=['F:\PL\BUPL\',area,'_roc_analysis\'];
    if ~exist(folderPrint,'dir')
        mkdir(folderPrint);
    end
    allChSessROCc=[];%compile ROC values across sessions and channels
    allChSessROCb=[];
    flipChs=[];
    excludeChs=[];
    markerCols=[0 1 1;1 0 1;0 1 0;0 0 0;1 0 0;0 0 1;0.5 0.5 1;1 0.5 0.5];
    sortCh=1;
    first10=0;
    if strcmp(area,'v1_3')
        sessionNums=127:130;
        channels=[7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];%Jack V1
        testContrasts=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
        excludeChs=[];
        flipChs=[];
    elseif strcmp(area,'v4_3')
        sessionNums=131:138;
        channels=[1:6 8 10 24 35 37 39 40 41 49 50 52:54 56];%Jack V4
        testContrasts=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
        excludeChs=[];%52 is mostly noise
        flipChs=[49];%49 has stimulus-evoked suppression
    end
    slpse_PLvsBUFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
    BUPLanalysisTypes=[{'BU'} {'PL'}];
    if readData==1
        for analysisTypeInd=1:2
            BUPLanalysisType=BUPLanalysisTypes{analysisTypeInd};
            sessROCvals=[];
            chROCvals=cell(length(channels),length(sessionNums));
            for j=1:length(sessionNums)
                rocvals=[];
                for condInd=1:length(testContrasts)
                    test_act=[];
                    sample_act=[];
                    for i=1:length(channels)%combine activity across channels                        
                        ch_test_act=[];
                        ch_sample_act=[];
                        matName=['F:\PL\BUPL\j_mat_files\',BUPLanalysisType,'\all_trials\',num2str(channels(i)),'_',num2str(sessionNums(j)),'.mat'];
                        if exist(matName,'file')
                            loadText=['load F:\PL\BUPL\j_mat_files\',BUPLanalysisType,'\all_trials\',num2str(channels(i)),'_',num2str(sessionNums(j)),'.mat'];
                            eval(loadText);
                            if size(matarray{condInd,2},1)~=size(matarray{condInd,4},1)
                                notEqual=1;
                                if size(matarray{condInd,2},1)<size(matarray{condInd,4},1)
                                    while size(matarray{condInd,2},1)~=size(matarray{condInd,4},1)
                                        matarray{condInd,2}=[NaN;matarray{condInd,2}];
                                        matarray{condInd,2}{1}=[];
                                        saveText=['save F:\PL\BUPL\j_mat_files\',BUPLanalysisType,'\all_trials\',num2str(channels(i)),'_',num2str(sessionNums(j)),'.mat matarray'];
                                        eval(saveText);
                                    end
                                elseif size(matarray{condInd,2},1)>size(matarray{condInd,4},1)
                                    while size(matarray{condInd,2},1)~=size(matarray{condInd,4},1)
                                        matarray{condInd,4}=[NaN;matarray{condInd,4}];
                                        matarray{condInd,4}{1}=[];
                                        saveText=['save F:\PL\BUPL\j_mat_files\',BUPLanalysisType,'\all_trials\',num2str(channels(i)),'_',num2str(sessionNums(j)),'.mat matarray'];
                                        eval(saveText);
                                    end
                                end
                            end
                            for trialInd=1:size(matarray{condInd,2},1)
                                test_act=[test_act length(matarray{condInd,4}{trialInd})];
                                sample_act=[sample_act length(matarray{condInd,2}{trialInd})];
                                ch_test_act=[ch_test_act length(matarray{condInd,4}{trialInd})];
                                ch_sample_act=[ch_sample_act length(matarray{condInd,2}{trialInd})];
                            end
                            ch_roc=sum(ch_test_act>ch_sample_act)/(sum(ch_test_act>ch_sample_act)+sum(ch_test_act<ch_sample_act));%new trial-wise method
                            chROCvals{i,j}=[chROCvals{i,j} ch_roc];
                        end
                    end
                    roc=sum(test_act>sample_act)/(sum(test_act>sample_act)+sum(test_act<sample_act));%new trial-wise method
                    rocvals=[rocvals roc];
                end
                sessROCvals=[sessROCvals;rocvals];
            end
            allROCvals{analysisTypeInd}=sessROCvals;
            allChROCvals{analysisTypeInd}=chROCvals;
        end
        saveText=['save F:\PL\BUPL\j_mat_files\',area,'cumulative_ROCvals_all.mat allROCvals'];
        eval(saveText)
        
        for analysisTypeInd=1:2
            chROCvals=allChROCvals{analysisTypeInd};
            BUPLanalysisType=BUPLanalysisTypes{analysisTypeInd};
            ROCvals=[];
            for j=1:length(sessionNums)
                for i=1:length(channels)%combine individual channel activity into a single matrix
                    if ~isempty(chROCvals{i,j})
                        ROCvals=[ROCvals;channels(i) sessionNums(j) chROCvals{i,j}];
                    end
                end
            end
            saveText=['save F:\PL\BUPL\j_mat_files\',BUPLanalysisType,'\',area,'ROCvals_all.mat ROCvals'];
            eval(saveText)
        end
    end
    xvals=0:1:testContrasts(1,end)+10;
    sampleContrast=30;
    if plotIndivChs==1
        for j=1:length(sessionNums)
            PSE=[];slope=[];
            ROC_PLvsBUFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
            sessionNum=sessionNums(j);
            % PL=[];
            if first10==1
                %             load(['F:\jack\before_300113\j_',area,'_roc_analysis_c\ROCvals_c_first10.mat'])
            else
                %             load(['F:\jack\before_300113\j_',area,'_roc_analysis_c\ROCvals_c_all.mat'])
                loadText=['load F:\PL\BUPL\j_mat_files\BU\',area,'ROCvals_all.mat'];
                eval(loadText);
            end
            PL=ROCvals;
            ind=find(PL(:,2)==sessionNum);
            PL=PL(ind,:);
            [temp order]=sort(PL(:,1));
            PL=PL(order,:);
            % BU=[];
            if first10==1
                %             load(['F:\jack\before_300113\j_',area,'_roc_analysis_b\ROCvals_b_first10.mat'])
            else
                %             load(['F:\jack\before_300113\j_',area,'_roc_analysis_b\ROCvals_b_all.mat'])
                loadText=['load F:\PL\BUPL\j_mat_files\PL\',area,'ROCvals_all.mat'];
                eval(loadText);
            end
            BU=ROCvals;
            ind=find(BU(:,2)==sessionNum);
            BU=BU(ind,:);
            [temp order]=sort(BU(:,1));
            BU=BU(order,:);
            countChs=0;
            for i=1:size(PL,1)
                ROCc=PL(i,3:end);%ROC values
                subplot(5,ceil(size(PL,1)/5),i);
                if channels(i)==52&&j==2
                    pauseHere=1;
                end
                %contrast discrimination task:
                plot(testContrasts(1,:),PL(i,3:end),'Marker','o','Color','b','LineStyle','none','MarkerFaceColor','b','MarkerEdgeColor','b');hold on
                if isempty(find(PL(i,1)==excludeChs))
                    countChs=countChs+1;
                    allChSessROCc=[allChSessROCc;ROCc];
                    if ~isempty(find(PL(i,1)==flipChs))%stim-induced suppression
                        allChSessROCc=[allChSessROCc;1-ROCc];
                    end
                    if sum(ROCc(1:3))<sum(ROCc(end-2:end))
                        X0=[2 30 0.5 0.1];
                    elseif sum(ROCc(1:3))>sum(ROCc(end-2:end))
                        X0=[-2 30 0.5 0.1];
                    end
                    X=fminsearch(@fit_weibull,X0,[],testContrasts(1,:),ROCc,[],'mle',[0 0 0 0],[],[0 0 0 0],[]);
                    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
                    xvalsFine=testContrasts(1):0.01:testContrasts(end);
                    yvalsFine=1-X(4)-X(3).*exp(-(xvalsFine./X(2)).^X(1));
                    if max(yvals)<0.5%out of range
                        PSE(countChs,1)=NaN;
                    elseif min(yvals)>0.5
                        PSE(countChs,1)=NaN;
                    else
                        diffTemp=yvalsFine-0.5;
                        [tempVal columnInd]=min(abs(diffTemp));
                        PSE(countChs,1)=xvalsFine(columnInd);
                        %PSE(countChs,1)=X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1));
                    end
                    if calculateTangent==0
                        slope(countChs,1)=X(1);
                    elseif calculateTangent==1
                        slope(countChs,1)=X(1)*X(3)*exp(-(sampleContrast/X(2))^X(1))*sampleContrast^(X(1)-1)*(1/X(2))^X(1);
                    end
                    line(PSE(countChs,1),0:0.01:1,'Color','b','LineStyle',':');
                    plot(xvals,yvals,'b');hold on
                    %orientation discrimination task:
                    plot(testContrasts(1,:),BU(i,3:end),'Marker','o','Color','r','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','r');hold on
                    ROCb=BU(i,3:end);%ROC values
                    if ~isempty(find(BU(i,1)==flipChs))%stim-induced suppression
                        allChSessROCb=[allChSessROCb;1-ROCb];
                    elseif isempty(find(BU(i,1)==excludeChs))
                        allChSessROCb=[allChSessROCb;ROCb];
                    end
                    if sum(ROCb(1:3))<sum(ROCb(end-2:end))
                        X0=[1 30 0.5 0.1];
                    elseif sum(ROCb(1:3))>sum(ROCb(end-2:end))
                        X0=[-1 30 0.5 0.1];
                    end
                    X=fminsearch(@fit_weibull,X0,[],testContrasts(1,:),ROCb,[],'mle',[0 0 0 0],[],[0 0 0 0],[]);
                    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
                    xvalsFine=testContrasts(1):0.01:testContrasts(end);
                    yvalsFine=1-X(4)-X(3).*exp(-(xvalsFine./X(2)).^X(1));
                    if max(yvals)<0.5%out of range
                        PSE(countChs,2)=NaN;
                    elseif min(yvals)>0.5
                        PSE(countChs,2)=NaN;
                    else
                        diffTemp=yvalsFine-0.5;
                        [tempVal columnInd]=min(abs(diffTemp));
                        PSE(countChs,2)=xvalsFine(columnInd);
                        %PSE(countChs,2)=X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1));
                    end
                    if calculateTangent==0
                        slope(countChs,2)=X(1);
                    elseif calculateTangent==1
                        slope(countChs,2)=X(1)*X(3)*exp(-(sampleContrast/X(2))^X(1))*sampleContrast^(X(1)-1)*(1/X(2))^X(1);
                    end
                    line(PSE(countChs,2),0:0.01:1,'Color','r','LineStyle',':');
                    plot(xvals,yvals,'r');
                    line(sampleContrast,0:0.01:1);
                    ylim([min([ROCc ROCb]),max([ROCc ROCb])]);
                    if strcmp(area,'v1_3')
                        xlim([0 100]);
                    elseif strcmp(area,'v4_3')
                        xlim([0 70]);
                    end
                    %     subplottitle=num2str(allMeanPerf(i,1));
                    subplottitle=num2str(BU(i,1));
                    title(subplottitle);
                    %     set(gca,'FontSize',[6],'YLim',[0,1.01],'XLim',[0,testContrasts(1:end)+10],'YTickMode','manual');%'YTick',[0.1]
                end
            end
            figure(ROC_PLvsBUFig);
            filename1=[folderPrint,num2str(sessionNum),'_PLvsBU_ROCs'];
            printtext=sprintf('print -dpng %s',filename1);
            eval(printtext);
            printtext=sprintf('print -depsc %s',filename1);
            eval(printtext);
            if sortCh==1
                if j==1%do sorting based on first session
                    [temp indp]=sort(PSE(:,1));
                end
                if size(PSE,1)~=length(indp)
                    indpTemp=[];
                    for rowInd=1:length(indp)
                        if find(indp(rowInd)==PL(:,1))
                            indpTemp=[indpTemp;indp(rowInd)];
                        end
                    end
                    PSEsorted=PSE(indpTemp,1:2);
                else
                    PSEsorted=PSE(indp,1:2);
                end
                if j==1
                    [temp inds]=sort(slope(:,1));
                end
                if size(slope,1)~=length(indp)
                    indsTemp=[];
                    for rowInd=1:length(indp)
                        if find(inds(rowInd)==PL(:,1))
                            indsTemp=[indsTemp;inds(rowInd)];
                        end
                    end
                    slopesorted=slope(indsTemp,1:2);
                else
                    slopesorted=slope(inds,1:2);
                end
            end
            PSEsorted=PSE;
            slopesorted=slope;
        end
        figure(slpse_PLvsBUFig)
        %     subplot(2,3,[1 2]);
        subplot(2,2,1);
        plot(1:size(PL,1),PSEsorted(:,1),'Marker','o','Color',markerCols(j,:),'LineStyle','none','MarkerFaceColor',markerCols(j,:),'MarkerEdgeColor',markerCols(j,:));hold on
        plot(1:size(BU,1),PSEsorted(:,2),'Marker','o','Color',markerCols(j,:),'LineStyle','none','MarkerEdgeColor',markerCols(j,:));hold on
        if strcmp(area,'v1_3')
            if first10==1
%                 ylim([10 120]);
            else
%                 ylim([20 90]);
            end
        elseif strcmp(area,'v4_3')
            if first10==1
                %ylim([]);
            else
                ylim([0 70]);
            end
        end
        %     subplot(2,3,3);
        subplot(2,2,2);
        plot(1:size(PL,1),(PSEsorted(:,1)-PSEsorted(:,2))./PSEsorted(:,2),'Marker','s','Color',markerCols(j,:),'LineStyle','none','MarkerFaceColor',markerCols(j,:),'MarkerEdgeColor',markerCols(j,:));hold on
        line([0 size(PL,1)],[0 0],'LineStyle',':');
        if strcmp(area,'v1_3')
            if first10==1
%                 ylim([-0.2 1.5]);
            else
%                 ylim([-0.6 0.2]);
            end
        elseif strcmp(area,'v4_3')
            if first10==1
                %ylim([]);
            else
%                 ylim([-2 2]);
            end
        end
        %     subplot(2,3,[4 5]);
        subplot(2,2,3);
        plot(1:size(PL,1),slopesorted(:,1),'Marker','o','Color',markerCols(j,:),'LineStyle','none','MarkerFaceColor',markerCols(j,:),'MarkerEdgeColor',markerCols(j,:));hold on
        plot(1:size(BU,1),slopesorted(:,2),'Marker','o','Color',markerCols(j,:),'LineStyle','none','MarkerEdgeColor',markerCols(j,:));hold on
        if strcmp(area,'v1_3')
            if first10==1
%                 ylim([0 7]);
            else
%                 ylim([-1 10]);
            end
        elseif strcmp(area,'v4_3')
            if first10==1
                %     ylim([]);
            else
%                 ylim([-3 45]);
            end
        end
        %     subplot(2,3,6);
        subplot(2,2,4);
        plot(1:size(PL,1),(slopesorted(:,1)-slopesorted(:,2))./slopesorted(:,2),'Marker','s','Color',markerCols(j,:),'LineStyle','none','MarkerFaceColor',markerCols(j,:),'MarkerEdgeColor',markerCols(j,:));hold on
        if strcmp(area,'v1_3')
            if first10==1
%                 ylim([-2 10]);
            else
%                 ylim([-1 1.5]);
            end
        elseif strcmp(area,'v4_3')
            if first10==1
                %     ylim([]);
            else
%                 ylim([-3 50]);
            end
        end
        line([0 length(channels)],[0 0],'LineStyle',':');
        figure(slpse_PLvsBUFig);
        filename1=[folderPrint,num2str(sessionNums(1)),'to',num2str(sessionNums(end)),'_PLvsBU_slpse'];
        printtext=sprintf('print -dpng %s',filename1);
        eval(printtext);
        printtext=sprintf('print -depsc %s',filename1);
        eval(printtext);
    end 
    
    loadText=['load F:\PL\BUPL\j_mat_files\',area,'cumulative_ROCvals_all.mat allROCvals'];
    eval(loadText)
    meanROCc=mean(allROCvals{2});
    meanROCb=mean(allROCvals{1});
    stdROCc=std(allROCvals{2});%PL task (attend-RF, contrast)
    stdROCb=std(allROCvals{1});%bottom-up attention task (attend-away, orientation)
%     meanROCc=mean(allChSessROCc,1);%PL task (attend-RF, contrast)
%     meanROCb=mean(allChSessROCb,1)%bottom-up attention task (attend-away, orientation)
%     stdROCc=std(allChSessROCc,1);%PL task (attend-RF, contrast)
%     stdROCb=std(allChSessROCb,1)%bottom-up attention task (attend-away, orientation)
    if samePlot==0
        allMeanROCFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
    elseif samePlot==1
        if areaCount==1
            allMeanROCFig=figure('Color',[1,1,1],'Units','Normalized','Position',[0.12, 0.08, 0.8, 0.8]);
        else
            figure(allMeanROCFig)
        end
        subplot(1,2,areaCount);
    end
    plot(testContrasts,meanROCc,'Marker','o','Color','b','LineStyle','none','MarkerFaceColor','b','MarkerEdgeColor','b');hold on
    errorbar(testContrasts,meanROCc,stdROCc,'b','LineStyle','none');%errorbar(X,Y,E) plots Y versus X with symmetric error bars 2*E(i) long
    if sum(meanROCc(1:3))<sum(meanROCc(end-2:end))
        X0=[2 30 0.2 0.1];
    elseif sum(meanROCc(1:3))>sum(meanROCc(end-2:end))
        X0=[-2 30 0.2 0.1];
    end
    X=fminsearch(@fit_weibull,X0,[],testContrasts(1,:),meanROCc,[],'mle',[0 0 0 0],[],[0 0 0 0],[]);
    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
    xvalsFine=testContrasts(1):0.01:testContrasts(end);
    yvalsFine=1-X(4)-X(3).*exp(-(xvalsFine./X(2)).^X(1));
    if max(yvals)<0.5%out of range
        PSEall(areaCount,1)=NaN;
    elseif min(yvals)>0.5
        PSEall(areaCount,1)=NaN;
    else
        diffTemp=yvalsFine-0.5;
        [tempVal columnInd]=min(abs(diffTemp));
        PSEall(areaCount,1)=xvalsFine(columnInd);
        %PSEall(areaCount,1)=X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1));
    end
    if calculateTangent==0
        slopeall(areaCount,1)=X(1);
    elseif calculateTangent==1
        slopeall(areaCount,1)=X(1)*X(3)*exp(-(sampleContrast/X(2))^X(1))*sampleContrast^(X(1)-1)*(1/X(2))^X(1);
    end
    line(PSEall(areaCount,1),0:0.01:1,'Color','b','LineStyle','-');
    plot([PSEall(areaCount,1) PSEall(areaCount,1)],[0 1],'Color','b','LineStyle','-');
    plot(xvals,yvals,'b');hold on
    plot(testContrasts,meanROCb,'Marker','o','Color','r','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','r');hold on
    errorbar(testContrasts,meanROCb,stdROCb,'r','LineStyle','none');%errorbar(X,Y,E) plots Y versus X with symmetric error bars 2*E(i) long
    if sum(meanROCb(1:3))<sum(meanROCb(end-2:end))
        X0=[2 30 0.2 0.1];
    elseif sum(meanROCb(1:3))>sum(meanROCb(end-2:end))
        X0=[-2 30 0.2 0.1];
    end
    X=fminsearch(@fit_weibull,X0,[],testContrasts(1,:),meanROCb,[],'mle',[0 0 0 0],[],[0 0 0 0],[]);
    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
    xvalsFine=testContrasts(1):0.01:testContrasts(end);
    yvalsFine=1-X(4)-X(3).*exp(-(xvalsFine./X(2)).^X(1));
    if max(yvals)<0.5%out of range
        PSEall(areaCount,2)=NaN;
    elseif min(yvals)>0.5
        PSEall(areaCount,2)=NaN;
    else
        diffTemp=yvalsFine-0.5;
        [tempVal columnInd]=min(abs(diffTemp));
        PSEall(areaCount,2)=xvalsFine(columnInd);
        %PSEall(areaCount,2)=X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1));
    end
    if calculateTangent==0
        slopeall(areaCount,2)=X(1);
    elseif calculateTangent==1
        slopeall(areaCount,2)=X(1)*X(3)*exp(-(sampleContrast/X(2))^X(1))*sampleContrast^(X(1)-1)*(1/X(2))^X(1);
    end
    line(PSEall(areaCount,2),0:0.01:1,'Color','r','LineStyle','-');
    plot([PSEall(areaCount,2) PSEall(areaCount,2)],[0 1],'Color','r','LineStyle','-');
    plot(xvals,yvals,'r');hold on
    if samePlot~=1||samePlot==1&&areaCount==1
        xlabel('contrast (%)');
        ylabel('PROBMAT');
    end
        if samePlot~=1
            title('Distributions of ROC values during contrast and orientation tasks');
            title('Distributions of ROC values during attend-RF (blue) and attend-away (red)');
        else
            title(areaTexts{areaCount});
        end
    figure(allMeanROCFig)
    if strcmp(area,'v4_3')
        ylim([0.3 0.7]);
    end
    figure(allMeanROCFig);
    if samePlot==0
        filename1=[folderPrint,num2str(sessionNums(1)),'to',num2str(sessionNums(end)),'_allChSessMeanROC_BUPL_',area,'large'];
        printtext=sprintf('print -dpng %s',filename1);
        eval(printtext);
        printtext=sprintf('print -depsc %s',filename1);
        eval(printtext);
    end 
    [h1(areaCount),p1(areaCount),ci,stats]=ttest(meanROCb,meanROCc);
    cis1{(areaCount)}=ci;
    stats1{(areaCount)}=stats;
    allROCvalsArr=[allROCvals{1} allROCvals{2}];
    attArr=[allROCvals{1}*0+1 allROCvals{2}*0+2];
    sessArr=[];
    for rowInd=1:size(allROCvals{1},1)
        sessArr=[sessArr;zeros(1,size(allROCvals{1},2)*2)+rowInd];
    end
    condArr=[];
    for rowInd=1:size(allROCvals{1},1)
        condArr=[condArr;1:length(testContrasts) 1:length(testContrasts)];
    end
    allROCvalsArr=reshape(allROCvalsArr,1,size(allROCvalsArr,1)*size(allROCvalsArr,2));
    attArr=reshape(attArr,1,size(attArr,1)*size(attArr,2));%1: attend away; 2: attend-RF
    sessArr=reshape(sessArr,1,size(sessArr,1)*size(sessArr,2));
    condArr=reshape(condArr,1,size(condArr,1)*size(condArr,2));
    [p,t,stats]=anovan(allROCvalsArr,{attArr,sessArr,condArr});
    p2{areaCount}=p;
    t2{areaCount}=t;
    stats2{areaCount}=stats;
    figure
    [buplANOVA{areaCount},m,h]=multcompare(stats,'dimension',[1 2])
    figure
    [buplANOVA{areaCount},m,h]=multcompare(stats,'dimension',[1 3])
    Fs(areaCount,1)=t{2,6};
end
if samePlot==1
    figure(allMeanROCFig); 
    subplot(1,2,2)
    xlim([0 70])
    ylim([0.25 0.7])
    folderPrint2='F:\PL\BUPL\';
    filename1=[folderPrint2,'allChSessMeanROC_BUPL_V4_3_V1_3_large'];
    printtext=sprintf('print -dpng %s',filename1);
    eval(printtext);
    printtext=sprintf('print -depsc %s',filename1);
    eval(printtext);
end
