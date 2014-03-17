function bj_roc_diff_figure(exampleFig,cutoff,animals,useISI,areas,excludeSuppressed,normalize)
%Written by Xing 01/10/13
%Generate figure for methods paper, comparing AUROC and COBAM methods.
%Read activity levels for all trials (pooled across sessions) for each
%channel, calculate individual channel AUROC and COBAM values.
%Rank channels in order of AUROC values, and then calculate cumulative
%AUROC/COBAM values when data from increasing numbers of channels included.

%Set useISI to 1: based on pre-test vs test, not on sample vs test.
%Set useISI to 0: sample vs test.

%Set exampleFig to 0 to plot ROC curves for new and old methods, set to 1
%to only plot example figures of distributions of stimulus-evoked activity
%and condition-dependent ROC curves.
%Set excludeSuppressed to 1 to exclude channels with stimulus-evoked
%suppression, i.e. blanco 13, 24, 42 and jack 49.
criterionType=0;
calculateTangent=1;
equalCount=0;
notEqualCount=0;
calcParams=0;
if useISI==1
    analysisType='ROC';
else
    analysisType='ROC_zero_one';
end
counter=0;%8 values- spiking B V4, B V1, J V4, J V1, MUA B V4, B V1, J V4, J V1
slopeNeuro_sglroc3=[];
PNE_sglroc3=[];
minRate_sglroc3=[];
maxRate_sglroc3=[];
diffPNE_sglroc3=[];
slopeNeuro=[];
PNE=[];
minRate=[];
maxRate=[];
diffPNE=[];
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
plotDiffC50_30=1;
if nargin<3||isempty(animals)
    animals=[{'blanco'} {'jack'}];
    % animals={'blanco'};
end
if nargin<4||isempty(areas)
    areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
    areas=[{'v4_1'} {'v1_1'}];
end
readData=0;
if readData==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            figROCconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(figROCconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
            [sampleContrasts testContrasts]=area_metadata(area);
            if strncmp(area,'v1_2',4)
                sampleContrasts=30;%
                testContrasts=testContrasts(2,:);
            end
            channels=main_channels(animal,area);
            sessionNums=main_raw_sessions_final(animal,area,[],0);
            if excludeSuppressed==1
                if strcmp(animal,'blanco')&&strcmp(area,'v4_1')%exclude channels 13, 24 and 42
                    channels=[1,2,3,4,7,12,14,18,20,22,33,34,36,37,38,40,49,50,51,52,53,54,55,57,59,60];
                elseif strcmp(animal,'jack')
                    channels=channels(channels~=49);
                end
            end
            for sampleContrastsInd=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(sampleContrastsInd);
                testContrast=testContrasts(sampleContrastsInd,:);
                %read in list of included channels
                if cutoff~=1
                    matname=['good_SNR_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];
                    pathname=fullfile(rootFolder,'PL','SNR',animal,matname);
                else
                    matname=['good_SNR_',area,'_',num2str(sampleContrast),'.mat'];
                    pathname=fullfile(rootFolder,'PL','SNR',animal,'cutoff_SNR_1',matname);
                end
                loadText=['load ',pathname,' includeSessionsAll'];
                eval(loadText);
                rocvals_sglroc3=zeros(1,length(testContrast));
                rocvals=zeros(1,length(testContrast));
                if useISI==0
                    allEpoch2=cell(length(channels),length(testContrast));%array to store activity for each ch and cond across sessions and trials
                elseif useISI==1
                    allEpoch3=cell(length(channels),length(testContrast));
                end
                allEpoch4=cell(length(channels),length(testContrast));
                for chInd=1:length(channels)
                    for i=1:length(sessionNums)%combine act across trials and across all sessions for each channel
                        for condInd=1:length(testContrast)
                            matPath=['F:\PL\sample_test_activity\',animal,'_',area,'\ch',num2str(channels(chInd)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat'];
                            loadText=['load ',matPath];
                            matExists=0;
                            if exist(matPath,'file')
                                matExists=1;
                            end
                            includeRows=includeSessionsAll(find(includeSessionsAll(:,1)==channels(chInd)),2);%include this session in analysis
                            includeRow=find(includeRows==sessionNums(i));
%                             if matExists==1&&~isempty(includeRow)&&includeCh
                                eval(loadText);
                                if useISI==0
                                    allEpoch2{chInd,condInd}=[allEpoch2{chInd,condInd} epoch2{condInd}];
                                elseif useISI==1
                                    allEpoch3{chInd,condInd}=[allEpoch3{chInd,condInd} epoch3{condInd}];%append act list across sessions
                                end
                                allEpoch4{chInd,condInd}=[allEpoch4{chInd,condInd} epoch4{condInd}];
                                if normalize==1
                                    maxAll=[];%find the highest firing rate across all conditions and trials, across both the sample and test presentation periods
                                    maxval=max(maxAll)/100;
                                else
                                    maxval=1;
                                end
%                             end
                            if chInd==1
                                numSessTrials(i,condInd)=size(epoch2{condInd},2);
                            end
                            if chInd==15
                                numSessTrials49(i,condInd)=size(epoch2{condInd},2);
                            end
                        end
                    end
                    for condInd=1:length(testContrast)
                        higherTestAct=0;
                        lowerTestAct=0;
                        for n=1:length(allEpoch4{chInd,condInd})%calculate COBAM for channel
                            if useISI==0                         
                                if criterionType==0%exclude trials where responses to sample and test are equal
                                    if allEpoch2{chInd,condInd}(n)<allEpoch4{chInd,condInd}(n)
                                        higherTestAct=higherTestAct+1;
                                    elseif allEpoch2{chInd,condInd}(n)>allEpoch4{chInd,condInd}(n)
                                        lowerTestAct=lowerTestAct+1;
                                    end
                                elseif criterionType==1%exclude trials where responses to sample and test are equal
                                    if allEpoch2{chInd,condInd}(n)<=allEpoch4{chInd,condInd}(n)
                                        higherTestAct=higherTestAct+1;
                                    elseif allEpoch2{chInd,condInd}(n)>allEpoch4{chInd,condInd}(n)
                                        lowerTestAct=lowerTestAct+1;
                                    end
                                elseif criterionType==2%exclude trials where responses to sample and test are equal
                                    if allEpoch2{chInd,condInd}(n)<allEpoch4{chInd,condInd}(n)
                                        higherTestAct=higherTestAct+1;
                                    elseif allEpoch2{chInd,condInd}(n)>=allEpoch4{chInd,condInd}(n)
                                        lowerTestAct=lowerTestAct+1;
                                    end
                                end
                                if allEpoch2{chInd,condInd}(n)==allEpoch4{chInd,condInd}(n)
                                    equalCount=equalCount+1;
                                else
                                    notEqualCount=notEqualCount+1;
                                end
                            elseif useISI==1
                                if allEpoch3{chInd,condInd}(n)<allEpoch4{chInd,condInd}(n)
                                    higherTestAct=higherTestAct+1;
                                elseif allEpoch3{chInd,condInd}(n)>allEpoch4{chInd,condInd}(n)
                                    lowerTestAct=lowerTestAct+1;
                                end
                            end
                        end
                        rocvals_sglroc3(chInd,condInd)=sglroc3(allEpoch4{chInd,condInd},allEpoch2{chInd,condInd});%old method
                        rocvals(chInd,condInd)=higherTestAct/(higherTestAct+lowerTestAct);%new method
                    end
                end
                allROCvals_sglroc3{areaInd,animalInd}=rocvals_sglroc3;
                allROCvals{areaInd,animalInd}=rocvals;
                allSubjectsEpoch4{areaInd,animalInd}=allEpoch4;
                allSubjectsEpoch2{areaInd,animalInd}=allEpoch2;
                allChannels{areaInd,animalInd}=channels;
                %             subplot(2,2,animalInd+2*(areaInd-1));
                for condInd=1:length(testContrast)
                    subplot(3,5,condInd);
                    plot(1:size(rocvals,1),rocvals(:,condInd),'bo');hold on
                    plot(1:size(rocvals,1),rocvals_sglroc3(:,condInd),'ro');
                    plot([0 size(rocvals,1)+1],[0.5 0.5],'k--');
                    xlim([0 size(rocvals,1)+1]);
                    title(num2str(condInd));
                end
            end
            if useISI==0
                subFolder='new_vs_old_individualchannels';
            elseif useISI==1
                subFolder='new_vs_old_useISI_individualchannels';
            end
            if excludeSuppressed
                subFolder=[subFolder,'_excludeSuppressed'];
            end
            if normalize
                subFolder=[subFolder,'_normalised'];
            end
            if criterionType==1
                subFolder=[subFolder,'_equalOrLower'];
            elseif criterionType==2
                subFolder=[subFolder,'_equalOrHigher'];
            end
            if useISI==0
                imagename=['all_ind_chs_ROCs_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
            elseif useISI==1
                imagename=['all_ind_chs_ROCs_old_new_useISI_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
            end
            folderpathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder);
            if ~exist(folderpathname,'dir')
                mkdir(folderpathname);
            end
            pathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder,imagename);
            printtext=sprintf('print -dpng %s.png',pathname);
            set(gcf,'PaperPositionMode','auto')
            eval(printtext);
        end
    end
    if excludeSuppressed==0
        saveText=['save ',rootFolder,'\PL\',analysisType,'\allIndividualChs_AUROC_COBAM_vals allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
        eval(saveText)
    elseif excludeSuppressed==1
        saveText=['save ',rootFolder,'\PL\',analysisType,'\allIndividualChs_AUROC_COBAM_vals_excludeSuppressed allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
        eval(saveText)
    end
end

%generate example figure from Jack V1_1 test contrast of 20%
if excludeSuppressed==0
    loadText=['load F:\PL\ROC_zero_one\allIndividualChs_AUROC_COBAM_vals allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
    eval(loadText)
elseif excludeSuppressed==1
    loadText=['load F:\PL\ROC_zero_one\allIndividualChs_AUROC_COBAM_vals_excludeSuppressed allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
    eval(loadText)
end
condInd=4;
allEpoch4=allSubjectsEpoch4{2,2};
allEpoch2=allSubjectsEpoch2{2,2};
[sortedROC_sglroc3 ind]=sort(1-allROCvals_sglroc3{2,2}(:,condInd));%lowest to highest
%[sortedROC_sglroc3ind]=sort(1-allROCvals_sglroc3{2,2}(:,condInd),1,'descend');%highest to lowest
channels=allChannels{2,2}(ind);
sampleContrasts=30;
cumulativeChAct2=zeros(1,length(allEpoch2{1,condInd}));
cumulativeChAct4=zeros(1,length(allEpoch2{1,condInd}));
for chInd=1:length(channels)%append activity from increasingly more channels
    cumulativeChAct2=cumulativeChAct2+allEpoch2{ind(chInd),condInd};
    cumulativeChAct4=cumulativeChAct4+allEpoch4{ind(chInd),condInd};
    higherTestAct=0;
    lowerTestAct=0;
    for n=1:length(cumulativeChAct4)%calculate COBAM for channel
        if useISI==0
            if criterionType==0%exclude trials where responses to sample and test are equal
                if cumulativeChAct2(n)<cumulativeChAct4(n)
                    higherTestAct=higherTestAct+1;
                elseif cumulativeChAct2(n)>cumulativeChAct4(n)
                    lowerTestAct=lowerTestAct+1;
                end
            elseif criterionType==1%equal to or lower
                if cumulativeChAct2(n)<=cumulativeChAct4(n)
                    higherTestAct=higherTestAct+1;
                elseif cumulativeChAct2(n)>cumulativeChAct4(n)
                    lowerTestAct=lowerTestAct+1;
                end
            elseif criterionType==2%equal to or higher
                if cumulativeChAct2(n)<cumulativeChAct4(n)
                    higherTestAct=higherTestAct+1;
                elseif cumulativeChAct2(n)>=cumulativeChAct4(n)
                    lowerTestAct=lowerTestAct+1;
                end
            end
        end
    end
    rocvals_sglroc3_cumulative(chInd)=sglroc3(cumulativeChAct4,cumulativeChAct2);%old method
    rocvals_cumulative(chInd)=higherTestAct/(higherTestAct+lowerTestAct);%new method    
end
figure
plot(1:size(allROCvals{2,2},1),allROCvals{2,2}(ind,4),'bo');hold on
plot(1:size(allROCvals{2,2},1),allROCvals_sglroc3{2,2}(ind,4),'ro');
plot(1:size(allROCvals{2,2},1),rocvals_cumulative,'bo','MarkerFaceColor','b');hold on
plot(1:size(allROCvals{2,2},1),rocvals_sglroc3_cumulative,'ro','MarkerFaceColor','r');
figure
plot(1:size(allROCvals{2,2},1),1-allROCvals{2,2}(ind,4),'MarkerEdgeColor',[0.6 0.6 0.6],'Marker','o','LineStyle','none','MarkerFaceColor',[0.6 0.6 0.6]);hold on
plot(1:size(allROCvals{2,2},1),1-allROCvals_sglroc3{2,2}(ind,4),'MarkerEdgeColor',[0.6 0.6 0.6],'Marker','o','LineStyle','none');
plot(1:size(allROCvals{2,2},1),1-rocvals_cumulative,'ko','MarkerFaceColor','k');hold on
plot(1:size(allROCvals{2,2},1),1-rocvals_sglroc3_cumulative,'ko');
ylabel('spiking AUROC or COBAM','FontAngle','italic');
xlabel('channel number','FontAngle','italic');

jackExample=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.45, 0.3]);
subplot(1,2,1);
plot(1:size(allROCvals{2,2},1),allROCvals{2,2}(ind,4),'MarkerEdgeColor',[0.6 0.6 0.6],'Marker','o','LineStyle','none','MarkerFaceColor',[0.6 0.6 0.6]);hold on
plot(1:size(allROCvals{2,2},1),allROCvals_sglroc3{2,2}(ind,4),'MarkerEdgeColor',[0.6 0.6 0.6],'Marker','o','LineStyle','none');
plot(1:size(allROCvals{2,2},1),rocvals_cumulative,'ko','MarkerFaceColor','k');hold on
plot(1:size(allROCvals{2,2},1),rocvals_sglroc3_cumulative,'ko');
line([0 25],[0.4 0.4],'Color',[0.6 0.6 0.6]);
ylabel('AUROC or COBAM','FontAngle','italic');
xlabel('number of pooled channels','FontAngle','italic');
title('Spiking activity');
set(gca, 'box', 'off')

%plots of range of individual channel data, plus cumulative values
%across all channels
figContrasts=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.25, 0.8]); %
set(figContrasts, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape');
if excludeSuppressed==0
    loadText=['load F:\PL\ROC_zero_one\allIndividualChs_AUROC_COBAM_vals allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
    eval(loadText)
elseif excludeSuppressed==1
    loadText=['load F:\PL\ROC_zero_one\allIndividualChs_AUROC_COBAM_vals_excludeSuppressed allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
    eval(loadText)
end
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        subplot(4,2,animalInd+2*(areaInd-1));
        maxCh_sglroc3=[];
        maxCh=[];
        allEpoch4=allSubjectsEpoch4{areaInd,animalInd};
        allEpoch2=allSubjectsEpoch2{areaInd,animalInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        if strncmp(area,'v1_2',4)
            sampleContrasts=30;%
            testContrasts=testContrasts(2,:);
        end
        for sampleContrastsInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastsInd);
            testContrast=testContrasts(sampleContrastsInd,:);
            colmapText=colormap(jet(size(testContrast,2)));
            colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
            plot([30 30],[0 1],'Color',[0.5 0.5 0.5]);
            for condInd=1:length(testContrast)
                cumulativeChAct2=zeros(1,length(allEpoch2{1,condInd}));
                cumulativeChAct4=zeros(1,length(allEpoch2{1,condInd}));
                if condInd<=7
                    maxCh_sglroc3(condInd)=min(allROCvals_sglroc3{areaInd,animalInd}(:,condInd));%find value of highest AUROC across channels
                    maxCh(condInd)=min(allROCvals{areaInd,animalInd}(:,condInd));%find value of highest COBAM across channels
                else
                    maxCh_sglroc3(condInd)=max(allROCvals_sglroc3{areaInd,animalInd}(:,condInd));%find value of highest AUROC across channels
                    maxCh(condInd)=max(allROCvals{areaInd,animalInd}(:,condInd));%find value of highest COBAM across channels
                end
                %plot(condInd,allROCvals{areaInd,animalInd}(:,condInd),'Color',colmapText(condInd,:),'Marker','x','MarkerSize',5);hold on
                plot(testContrast(condInd),allROCvals{areaInd,animalInd}(:,condInd),'Color',colmapText(condInd,:),'Marker','.','MarkerSize',3);hold on
                plot(testContrast(condInd),allROCvals{areaInd,animalInd}(:,condInd),'b.','MarkerSize',3);hold on
                %plot(testContrast(condInd),mean(allROCvals{areaInd,animalInd}(:,condInd)),'Color',colmapText(condInd,:),'Marker','x');hold on
                %plot(condInd+0.4,allROCvals_sglroc3{areaInd,animalInd}(:,condInd),'Color',colmapText(condInd,:),'Marker','.');hold on
                iqr_roc(condInd)=iqr(allROCvals{areaInd,animalInd}(:,condInd));
                iqr_roc_upper(condInd)=mean(allROCvals{areaInd,animalInd}(:,condInd))+iqr_roc(condInd);
                iqr_roc_lower(condInd)=mean(allROCvals{areaInd,animalInd}(:,condInd))-iqr_roc(condInd);
                iqr_roc_sglroc3(condInd)=iqr(allROCvals_sglroc3{areaInd,animalInd}(:,condInd));
                iqr_roc_sglroc3_upper(condInd)=mean(allROCvals_sglroc3{areaInd,animalInd}(:,condInd))+iqr_roc_sglroc3(condInd);
                iqr_roc_sglroc3_lower(condInd)=mean(allROCvals_sglroc3{areaInd,animalInd}(:,condInd))-iqr_roc_sglroc3(condInd);
                for chInd=1:size(allEpoch2,1)%read activity across cells corresponding to channels
                    cumulativeChAct2=cumulativeChAct2+allEpoch2{chInd,condInd};
                    cumulativeChAct4=cumulativeChAct4+allEpoch4{chInd,condInd};
                end
                higherTestAct=0;
                lowerTestAct=0;
                for n=1:length(cumulativeChAct4)%calculate COBAM for channel
                    if useISI==0
                        if criterionType==0%exclude trials where responses to sample and test are equal
                            if cumulativeChAct2(n)<cumulativeChAct4(n)
                                higherTestAct=higherTestAct+1;
                            elseif cumulativeChAct2(n)>cumulativeChAct4(n)
                                lowerTestAct=lowerTestAct+1;
                            end
                        elseif criterionType==1%equal to or lower
                            if cumulativeChAct2(n)<=cumulativeChAct4(n)
                                higherTestAct=higherTestAct+1;
                            elseif cumulativeChAct2(n)>cumulativeChAct4(n)
                                lowerTestAct=lowerTestAct+1;
                            end
                        elseif criterionType==2%equal to or higher
                            if cumulativeChAct2(n)<cumulativeChAct4(n)
                                higherTestAct=higherTestAct+1;
                            elseif cumulativeChAct2(n)>=cumulativeChAct4(n)
                                lowerTestAct=lowerTestAct+1;
                            end
                        end
                    end
                end
                rocvals_sglroc3_allCumulative(condInd)=sglroc3(cumulativeChAct4,cumulativeChAct2);%old method
                rocvals_allCumulative(condInd)=higherTestAct/(higherTestAct+lowerTestAct);%new method
                plot(testContrast(condInd),rocvals_sglroc3_allCumulative(condInd),'ro','MarkerSize',3);
                plot(testContrast(condInd),rocvals_allCumulative(condInd),'bo','MarkerSize',3);
            end
            fillCoords=[iqr_roc_upper fliplr(iqr_roc_lower)];
            fillCoords_sglroc3=[iqr_roc_sglroc3_upper fliplr(iqr_roc_sglroc3_lower)];
            xfillCoords=[testContrast fliplr(testContrast)];
            patch(xfillCoords,fillCoords,'b','EdgeColor','none');
            patch(xfillCoords,fillCoords_sglroc3,'r','EdgeColor','none');alpha(0.3)
            xlimVals=get(gca,'XLim');
%             plot([xlimVals(1) xlimVals(2)],[0.5 0.5],'k--');
            if areaInd==2
                xlabel('test contrast (%)','FontAngle','italic');
            end
            if animalInd==1
                ylabel('AUROC or COBAM value','FontAngle','italic');
            end
            if calcParams==1
                [slopeNeuro_sglroc3,PNE_sglroc3,diffPNE_sglroc3,minRate_sglroc3,maxRate_sglroc3,chSSE_sglroc3]=weibull_fitting(rocvals_sglroc3_allCumulative,sampleContrast,testContrast,'old',counter,slopeNeuro_sglroc3,[],PNE_sglroc3,minRate_sglroc3,maxRate_sglroc3,diffPNE_sglroc3,0,calculateTangent,0);
                [slopeNeuro,PNE,diffPNE,minRate,maxRate,chSSE]=weibull_fitting(rocvals_allCumulative,sampleContrast,testContrast,'old',counter,slopeNeuro,[],PNE,minRate,maxRate,diffPNE,0,calculateTangent,0);
            end
            xlim([testContrast(1) testContrast(end)]);
            ylim([0 1]);
            axis square
            set(gca, 'box', 'off')
            counter=counter+1;
        end
    end
end
equalCount
notEqualCount
%plot2svg%xport as an SVG file and open in Illustrator, then save as an EPS/TIF

%defunct: reads data in which identity of trials is scrambled
convertMUAfiles=0;
if convertMUAfiles==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            channels=main_channels(animal,area);
            sessionNums=main_raw_sessions_final(animal,area,[],0);
            if excludeSuppressed==1
                if strcmp(animal,'blanco')&&strcmp(area,'v4_1')%exclude channels 13, 24 and 42
                    channels=[1,2,3,4,7,12,14,18,20,22,33,34,36,37,38,40,49,50,51,52,53,54,55,57,59,60];
                elseif strcmp(animal,'jack')
                    channels=channels(channels~=49);
                end
            end
            for sampleContrastsInd=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(sampleContrastsInd);
                testContrast=testContrasts(sampleContrastsInd,:);
                %read in list of included channels
                if useISI==0
                    allEpoch2=cell(length(channels),length(testContrast));%array to store activity for each ch and cond across sessions and trials
                elseif useISI==1
                    allEpoch3=cell(length(channels),length(testContrast));
                end
                allEpoch4=cell(length(channels),length(testContrast));
                for chInd=1:length(channels)
                    for i=1:length(sessionNums)%combine act across trials and across all sessions for each channel
                        loadText=['load I:\MUA_files_AUROC_COBAM\',animal,'\',num2str(sessionNums(i)),'\FR_MUA_50_533\FR_MUA_',num2str(channels(chInd)),'_50_533'];
                        eval(loadText)
                        for condInd=1:length(testContrast)
                            epoch2{condInd,1}=[FR_SAMPLE_CORR(1:CND_Trial_Counter_CORR(condInd)) FR_SAMPLE_ERR(1:CND_Trial_Counter_ERR(condInd))];
                            epoch4{condInd,1}=[FR_TEST_CORR(1:CND_Trial_Counter_CORR(condInd)) FR_TEST_ERR(1:CND_Trial_Counter_ERR(condInd))];
                        end
                        saveText=['save F:\PL\sample_test_MUA\',animal,'_',area,'\ch',num2str(channels(chInd)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat epoch2 epoch4'];
                        eval(saveText)
                    end
                end
            end
        end
    end
end

convertMUAfiles2=0;
if convertMUAfiles2==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            channels=main_channels(animal,area);
            sessionNums=main_raw_sessions_final(animal,area,[],0);
            MUAsessionNumsV4{1}=[304:308 311 313 314 316:318 320:324 327:359];%included blanco MUA sessions
            includeSessionsInd=[4:8 10:13 19:31];
            if excludeSuppressed==1
                if strcmp(animal,'blanco')&&strcmp(area,'v4_1')%exclude channels 13, 24 and 42
                    channels=[1,2,3,4,7,12,14,18,20,22,33,34,36,37,38,40,49,50,51,52,53,54,55,57,59,60];
                end
                if strcmp(animal,'jack')
                    channels=channels(channels~=49);
                end
            end
            for sampleContrastsInd=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(sampleContrastsInd);
                testContrast=testContrasts(sampleContrastsInd,:);
                %read in list of included channels
                if useISI==0
                    allEpoch2=cell(length(channels),length(testContrast));%array to store activity for each ch and cond across sessions and trials
                elseif useISI==1
                    allEpoch3=cell(length(channels),length(testContrast));
                end
                allEpoch4=cell(length(channels),length(testContrast));
                loadText=['load I:\MUA_files_AUROC_COBAM\',animal,'_',area,'_populationRes'];
                eval(loadText)
                excludeTrialInd=cell(length(testContrast),length(sessionNums));
                numTrialsTotal=cell(length(testContrast),length(sessionNums));
                for chInd=1:length(channels)
                    for i=1:length(sessionNums)%combine act across trials and across all sessions for each channel
                        for condInd=1:length(testContrast)
                            if strcmp(area,'v4_1')&&strcmp(animal,'blanco')
                                epoch2{condInd,1}=[FR_SAMPLE_CORR{condInd,chInd,includeSessionsInd(i)} FR_SAMPLE_ERR{condInd,chInd,includeSessionsInd(i)}];%for V4, add three to the value of i, because matrix includes three sessions at beginning with 304 and 306 (8 and 12 conditions) and 306 (excluded), and also 342 at end w horizontal Gabor
                                epoch4{condInd,1}=[FR_TEST_CORR{condInd,chInd,includeSessionsInd(i)} FR_TEST_ERR{condInd,chInd,includeSessionsInd(i)}];
                            else
                                epoch2{condInd,1}=[FR_SAMPLE_CORR{condInd,chInd,i} FR_SAMPLE_ERR{condInd,chInd,i}];
                                epoch4{condInd,1}=[FR_TEST_CORR{condInd,chInd,i} FR_TEST_ERR{condInd,chInd,i}];
                            end
                            numTrialsTotal{condInd,i}=[numTrialsTotal{condInd,i} length(epoch2{condInd})];
                        end
                        saveText=['save F:\PL\sample_test_MUA\',animal,'_',area,'\ch',num2str(channels(chInd)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat epoch2 epoch4'];
                        eval(saveText)
                        for condInd=1:length(testContrast)%compile list of trials that get excluded for each session, across channels
                            excludeTrialInd{condInd,i}=[excludeTrialInd{condInd,i} find(isnan(epoch2{condInd})) find(isnan(epoch4{condInd}))];%combine exclusion list across sample and test presentations
                        end
                    end
                end
                for rowInd=1:length(testContrast)
                    for columnInd=1:length(sessionNums)
                        for chIndex=1:length(channels)-1
                            if numTrialsTotal{rowInd,columnInd}(chIndex+1)~=numTrialsTotal{rowInd,columnInd}(1)
                                numTrialsTotal{rowInd,columnInd}(chIndex+1)
                                numTrialsTotal{rowInd,columnInd}(1)
                            end
                        end
                    end
                end
                %look through arrays for NaN values and exclude trials from
                %all channels if at least 1 channel has a NaN entry for
                %that trial:
                excludeTrialFraction=cell(length(channels),length(sessionNums));
                for i=1:length(sessionNums)
                    for condInd=1:length(testContrast)
                        excludeTrialInd{condInd,i}=unique(excludeTrialInd{condInd,i});
                        for chInd=1:length(channels)
                            loadText=['load F:\PL\sample_test_MUA\',animal,'_',area,'\ch',num2str(channels(chInd)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat epoch2 epoch4'];
                            eval(loadText)
                            if chInd==1
                                excludeTrialLogicalInd=zeros(1,length(epoch2{condInd}));
                                for trialCount=1:length(epoch2{condInd})
                                    if find(trialCount==excludeTrialInd{condInd,i})
                                        excludeTrialLogicalInd(trialCount)=1;
                                    end
                                end
                                excludeTrialFraction{condInd,i}=length(excludeTrialInd{condInd,i})/length(epoch2{condInd});
                            end
                            if chInd==length(channels)
                                pauseHere=1;
                            end
                            if chInd==1
                                numTrials1=length(epoch2{condInd});
                            elseif chInd==length(channels)
                                if numTrials1~=length(epoch2{condInd});
                                    sessionNums(i)
                                    condInd
                                end
                            end
                            epoch2{condInd}=epoch2{condInd}(~excludeTrialLogicalInd);
                            epoch4{condInd}=epoch4{condInd}(~excludeTrialLogicalInd);
                            if chInd==1
                                numTrials2=length(epoch2{condInd});
                            elseif chInd==length(channels)
                                if numTrials2~=length(epoch2{condInd});
                                    sessionNums(i)
                                    condInd
                                end
                            end
                            saveText=['save F:\PL\sample_test_MUA\',animal,'_',area,'\ch',num2str(channels(chInd)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat epoch2 epoch4'];
                            eval(saveText)
                        end
                    end
                end
                saveText=['save F:\PL\sample_test_MUA\',animal,'_',area,'_excludeTrialFraction.mat excludeTrialFraction'];
                eval(saveText)
            end
        end
    end
end

readData2=0;
if readData2==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            figROCconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(figROCconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
            [sampleContrasts testContrasts]=area_metadata(area);
            if strncmp(area,'v1_2',4)
                sampleContrasts=30;%
                testContrasts=testContrasts(2,:);
            end
            channels=main_channels(animal,area);
            sessionNums=main_raw_sessions_final(animal,area,[],0);
            if excludeSuppressed==1
                if strcmp(animal,'blanco')&&strcmp(area,'v4_1')%exclude channels 13, 24 and 42
                    channels=[1,2,3,4,7,12,14,18,20,22,33,34,36,37,38,40,49,50,51,52,53,54,55,57,59,60];
                elseif strcmp(animal,'jack')
                    channels=channels(channels~=49);
                end
            end
            for sampleContrastsInd=1:length(sampleContrasts)
                sampleContrast=sampleContrasts(sampleContrastsInd);
                testContrast=testContrasts(sampleContrastsInd,:);
                %read in list of included channels
                rocvals_sglroc3=zeros(1,length(testContrast));
                rocvals=zeros(1,length(testContrast));
                if useISI==0
                    allEpoch2=cell(length(channels),length(testContrast));%array to store activity for each ch and cond across sessions and trials
                elseif useISI==1
                    allEpoch3=cell(length(channels),length(testContrast));
                end
                allEpoch4=cell(length(channels),length(testContrast));
                for chInd=1:length(channels)
                    for i=1:length(sessionNums)%combine act across trials and across all sessions for each channel
                        for condInd=1:length(testContrast)
                            matPath=['F:\PL\sample_test_MUA\',animal,'_',area,'\ch',num2str(channels(chInd)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat'];
                            loadText=['load ',matPath];
                            matExists=0;
                            if exist(matPath,'file')
                                matExists=1;
                            end
%                             if matExists==1&&~isempty(includeRow)&&includeCh
                                eval(loadText);
                                if useISI==0
                                    allEpoch2{chInd,condInd}=[allEpoch2{chInd,condInd} epoch2{condInd}];
                                elseif useISI==1
                                    allEpoch3{chInd,condInd}=[allEpoch3{chInd,condInd} epoch3{condInd}];%append act list across sessions
                                end
                                allEpoch4{chInd,condInd}=[allEpoch4{chInd,condInd} epoch4{condInd}];
                                if normalize==1
                                    maxAll=[];%find the highest firing rate across all conditions and trials, across both the sample and test presentation periods
                                    maxval=max(maxAll)/100;
                                else
                                    maxval=1;
                                end
%                             end
                            if chInd==1
                                numSessTrials(i,condInd)=size(epoch2{condInd},2);
                            end
                            if chInd==15
                                numSessTrials49(i,condInd)=size(epoch2{condInd},2);
                            end
                        end
                    end
                    for condInd=1:length(testContrast)
                        higherTestAct=0;
                        lowerTestAct=0;
                        for n=1:length(allEpoch4{chInd,condInd})%calculate COBAM for channel
                            if useISI==0
                                if criterionType==0%exclude trials where responses to sample and test are equal
                                    if allEpoch2{chInd,condInd}(n)<allEpoch4{chInd,condInd}(n)
                                        higherTestAct=higherTestAct+1;
                                    elseif allEpoch2{chInd,condInd}(n)>allEpoch4{chInd,condInd}(n)
                                        lowerTestAct=lowerTestAct+1;
                                    end
                                elseif criterionType==1%equal to or lower
                                    if allEpoch2{chInd,condInd}(n)<=allEpoch4{chInd,condInd}(n)
                                        higherTestAct=higherTestAct+1;
                                    elseif allEpoch2{chInd,condInd}(n)>allEpoch4{chInd,condInd}(n)
                                        lowerTestAct=lowerTestAct+1;
                                    end
                                elseif criterionType==2%equal to or higher
                                    if allEpoch2{chInd,condInd}(n)<allEpoch4{chInd,condInd}(n)
                                        higherTestAct=higherTestAct+1;
                                    elseif allEpoch2{chInd,condInd}(n)>=allEpoch4{chInd,condInd}(n)
                                        lowerTestAct=lowerTestAct+1;
                                    end
                                end
                                if allEpoch2{chInd,condInd}(n)==allEpoch4{chInd,condInd}(n)
                                    equalCount=equalCount+1;
                                else
                                    notEqualCount=notEqualCount+1;
                                end
                            elseif useISI==1
                                if allEpoch3{chInd,condInd}(n)<allEpoch4{chInd,condInd}(n)
                                    higherTestAct=higherTestAct+1;
                                elseif allEpoch3{chInd,condInd}(n)>allEpoch4{chInd,condInd}(n)
                                    lowerTestAct=lowerTestAct+1;
                                end
                            end
                        end
                        %traditional ROC method:
                        rocvals_sglroc3(chInd,condInd)=rocMUA(allEpoch4{chInd,condInd},allEpoch2{chInd,condInd});%traditional method, can handle negative MUA values
                        %rocvals_sglroc3(chInd,condInd)=sglroc3(allEpoch4{chInd,condInd},allEpoch2{chInd,condInd});%old method, unable to handle negative MUA values
                        rocvals(chInd,condInd)=higherTestAct/(higherTestAct+lowerTestAct);%new method
                    end
                end
                allROCvals_sglroc3{areaInd,animalInd}=rocvals_sglroc3;
                allROCvals{areaInd,animalInd}=rocvals;
                allSubjectsEpoch4{areaInd,animalInd}=allEpoch4;
                allSubjectsEpoch2{areaInd,animalInd}=allEpoch2;
                allChannels{areaInd,animalInd}=channels;
                %             subplot(2,2,animalInd+2*(areaInd-1));
                for condInd=1:length(testContrast)
                    subplot(3,5,condInd);
                    plot(1:size(rocvals,1),rocvals(:,condInd),'bo');hold on
                    plot(1:size(rocvals,1),rocvals_sglroc3(:,condInd),'ro');
                    plot([0 size(rocvals,1)+1],[0.5 0.5],'k--');
                    xlim([0 size(rocvals,1)+1]);
                    title(num2str(condInd));
                end
            end
            if useISI==0
                subFolder='new_vs_old_individualchannels_MUA';
            elseif useISI==1
                subFolder='new_vs_old_useISI_individualchannels_MUA';
            end
            if excludeSuppressed
                subFolder=[subFolder,'_excludeSuppressed'];
            end
            if normalize
                subFolder=[subFolder,'_normalised'];
            end
            if criterionType==1
                subFolder=[subFolder,'_equalOrLower'];
            elseif criterionType==2
                subFolder=[subFolder,'_equalOrHigher'];
            end
            if useISI==0
                imagename=['all_ind_chs_ROCs_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
            elseif useISI==1
                imagename=['all_ind_chs_ROCs_old_new_useISI_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
            end
            folderpathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder);
            if ~exist(folderpathname,'dir')
                mkdir(folderpathname);
            end
            pathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder,imagename);
            printtext=sprintf('print -dpng %s.png',pathname);
            set(gcf,'PaperPositionMode','auto')
            eval(printtext);
        end
    end
    if excludeSuppressed==0
        saveText=['save ',rootFolder,'\PL\',analysisType,'\allIndividualChs_AUROC_COBAM_vals_MUA allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
        eval(saveText)
    elseif excludeSuppressed==1
        saveText=['save ',rootFolder,'\PL\',analysisType,'\allIndividualChs_AUROC_COBAM_vals_excludeSuppressed_MUA allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
        eval(saveText)
    end
end

%generate MUA example figure from Jack V1_1 test contrast of 20%
if excludeSuppressed==0
    loadText=['load F:\PL\ROC_zero_one\allIndividualChs_AUROC_COBAM_vals_MUA allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
    eval(loadText)
elseif excludeSuppressed==1
    loadText=['load F:\PL\ROC_zero_one\allIndividualChs_AUROC_COBAM_vals_excludeSuppressed_MUA allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
    eval(loadText)
end
condInd=4;
allEpoch4=allSubjectsEpoch4{2,2};
allEpoch2=allSubjectsEpoch2{2,2};
[sortedROC_sglroc3 ind]=sort(1-allROCvals_sglroc3{2,2}(:,condInd));%lowest to highest
%[sortedROC_sglroc3ind]=sort(1-allROCvals_sglroc3{2,2}(:,condInd),1,'descend');%highest to lowest
channels=allChannels{2,2}(ind);
sampleContrasts=30;
rocvals_sglroc3_cumulative=[];
rocvals_cumulative=[];
cumulativeChAct2=zeros(1,length(allEpoch2{1,condInd}));
cumulativeChAct4=zeros(1,length(allEpoch2{1,condInd}));
for chInd=1:length(channels)%append activity from increasingly more channels
    cumulativeChAct2=cumulativeChAct2+allEpoch2{ind(chInd),condInd};
    cumulativeChAct4=cumulativeChAct4+allEpoch4{ind(chInd),condInd};
    higherTestAct=0;
    lowerTestAct=0;
    for n=1:length(cumulativeChAct4)%calculate COBAM for channel
        if useISI==0
            if criterionType==0%exclude trials where responses to sample and test are equal
                if cumulativeChAct2(n)<cumulativeChAct4(n)
                    higherTestAct=higherTestAct+1;
                elseif cumulativeChAct2(n)>cumulativeChAct4(n)
                    lowerTestAct=lowerTestAct+1;
                end
            elseif criterionType==1%equal to or lower
                if cumulativeChAct2(n)<=cumulativeChAct4(n)
                    higherTestAct=higherTestAct+1;
                elseif cumulativeChAct2(n)>cumulativeChAct4(n)
                    lowerTestAct=lowerTestAct+1;
                end
            elseif criterionType==2%equal to or higher
                if cumulativeChAct2(n)<cumulativeChAct4(n)
                    higherTestAct=higherTestAct+1;
                elseif cumulativeChAct2(n)>=cumulativeChAct4(n)
                    lowerTestAct=lowerTestAct+1;
                end
            end
        end
    end
    rocvals_sglroc3_cumulative(chInd)=rocMUA(cumulativeChAct4,cumulativeChAct2);%traditional method, can handle negative MUA values
    %rocvals_sglroc3_cumulative(chInd)=sglroc3(cumulativeChAct4,cumulativeChAct2);%old method
    rocvals_cumulative(chInd)=higherTestAct/(higherTestAct+lowerTestAct);%new method    
end
figure
plot(1:size(allROCvals{2,2},1),allROCvals{2,2}(ind,4),'bo');hold on
plot(1:size(allROCvals{2,2},1),allROCvals_sglroc3{2,2}(ind,4),'ro');
plot(1:size(allROCvals{2,2},1),rocvals_cumulative,'bo','MarkerFaceColor','b');hold on
plot(1:size(allROCvals{2,2},1),rocvals_sglroc3_cumulative,'ro','MarkerFaceColor','r');
figure
plot(1:size(allROCvals{2,2},1),1-allROCvals{2,2}(ind,4),'MarkerEdgeColor',[0.6 0.6 0.6],'Marker','o','LineStyle','none','MarkerFaceColor',[0.6 0.6 0.6]);hold on
plot(1:size(allROCvals{2,2},1),1-allROCvals_sglroc3{2,2}(ind,4),'MarkerEdgeColor',[0.6 0.6 0.6],'Marker','o','LineStyle','none');
plot(1:size(allROCvals{2,2},1),1-rocvals_cumulative,'ko','MarkerFaceColor','k');hold on
plot(1:size(allROCvals{2,2},1),1-rocvals_sglroc3_cumulative,'ko');
ylabel('MUA AUROC or COBAM','FontAngle','italic');
xlabel('channel number','FontAngle','italic');

figure(jackExample);
subplot(1,2,2);
plot(1:size(allROCvals{2,2},1),allROCvals{2,2}(ind,4),'MarkerEdgeColor',[0.6 0.6 0.6],'Marker','o','LineStyle','none','MarkerFaceColor',[0.6 0.6 0.6]);hold on
plot(1:size(allROCvals{2,2},1),allROCvals_sglroc3{2,2}(ind,4),'MarkerEdgeColor',[0.6 0.6 0.6],'Marker','o','LineStyle','none');
plot(1:size(allROCvals{2,2},1),rocvals_cumulative,'ko','MarkerFaceColor','k');hold on
plot(1:size(allROCvals{2,2},1),rocvals_sglroc3_cumulative,'ko');
line([0 25],[0.4 0.4],'Color',[0.6 0.6 0.6]);
ylim([0 0.4]);
xlabel('number of pooled channels','FontAngle','italic');
title('MUA');
set(gca, 'box', 'off')

%plots of range of individual channel MUA data, plus cumulative values
%across all channels
figure(figContrasts);
if excludeSuppressed==0
    loadText=['load F:\PL\ROC_zero_one\allIndividualChs_AUROC_COBAM_vals_MUA allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
    eval(loadText)
elseif excludeSuppressed==1
    loadText=['load F:\PL\ROC_zero_one\allIndividualChs_AUROC_COBAM_vals_excludeSuppressed_MUA allROCvals_sglroc3 allROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels'];
    eval(loadText)
end
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        subplot(4,2,animalInd+2*(areaInd-1)+4);
        maxCh_sglroc3=[];
        maxCh=[];
        allEpoch2=allSubjectsEpoch2{areaInd,animalInd};
        allEpoch4=allSubjectsEpoch4{areaInd,animalInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        if strncmp(area,'v1_2',4)
            sampleContrasts=30;%
            testContrasts=testContrasts(2,:);
        end
        for sampleContrastsInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastsInd);
            testContrast=testContrasts(sampleContrastsInd,:);
            colmapText=colormap(jet(size(testContrast,2)));
            colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
            plot([30 30],[0 1],'Color',[0.5 0.5 0.5]);
            for condInd=1:length(testContrast)
                cumulativeChAct2=zeros(1,length(allEpoch2{1,condInd}));
                cumulativeChAct4=zeros(1,length(allEpoch4{1,condInd}));
                if condInd<=7
                    maxCh_sglroc3(condInd)=min(allROCvals_sglroc3{areaInd,animalInd}(:,condInd));%find value of highest AUROC across channels
                    maxCh(condInd)=min(allROCvals{areaInd,animalInd}(:,condInd));%find value of highest COBAM across channels
                else
                    maxCh_sglroc3(condInd)=max(allROCvals_sglroc3{areaInd,animalInd}(:,condInd));%find value of highest AUROC across channels
                    maxCh(condInd)=max(allROCvals{areaInd,animalInd}(:,condInd));%find value of highest COBAM across channels
                end
                %plot(condInd,allROCvals{areaInd,animalInd}(:,condInd),'Color',colmapText(condInd,:),'Marker','x','MarkerSize',5);hold on
                plot(testContrast(condInd),allROCvals{areaInd,animalInd}(:,condInd),'Color',colmapText(condInd,:),'Marker','.','MarkerSize',3);hold on
                plot(testContrast(condInd),allROCvals{areaInd,animalInd}(:,condInd),'b.','MarkerSize',3);hold on
                %plot(testContrast(condInd),mean(allROCvals{areaInd,animalInd}(:,condInd)),'Color',colmapText(condInd,:),'Marker','x');hold on
                %plot(condInd+0.4,allROCvals_sglroc3{areaInd,animalInd}(:,condInd),'Color',colmapText(condInd,:),'Marker','.');hold on
                iqr_roc(condInd)=iqr(allROCvals{areaInd,animalInd}(:,condInd));
                iqr_roc_upper(condInd)=mean(allROCvals{areaInd,animalInd}(:,condInd))+iqr_roc(condInd);
                iqr_roc_lower(condInd)=mean(allROCvals{areaInd,animalInd}(:,condInd))-iqr_roc(condInd);
                iqr_roc_sglroc3(condInd)=iqr(allROCvals_sglroc3{areaInd,animalInd}(:,condInd));
                iqr_roc_sglroc3_upper(condInd)=mean(allROCvals_sglroc3{areaInd,animalInd}(:,condInd))+iqr_roc_sglroc3(condInd);
                iqr_roc_sglroc3_lower(condInd)=mean(allROCvals_sglroc3{areaInd,animalInd}(:,condInd))-iqr_roc_sglroc3(condInd);
                for chInd=1:size(allEpoch2,1)%read activity across cells corresponding to channels
                    cumulativeChAct2=cumulativeChAct2+allEpoch2{chInd,condInd};
                    cumulativeChAct4=cumulativeChAct4+allEpoch4{chInd,condInd};
                end
                higherTestAct=0;
                lowerTestAct=0;
                for n=1:length(cumulativeChAct4)%calculate COBAM for channel
                    if useISI==0
                        if criterionType==0%exclude trials where responses to sample and test are equal
                            if cumulativeChAct2(n)<cumulativeChAct4(n)
                                higherTestAct=higherTestAct+1;
                            elseif cumulativeChAct2(n)>cumulativeChAct4(n)
                                lowerTestAct=lowerTestAct+1;
                            end
                        elseif criterionType==1%equal to or lower
                            if cumulativeChAct2(n)<=cumulativeChAct4(n)
                                higherTestAct=higherTestAct+1;
                            elseif cumulativeChAct2(n)>cumulativeChAct4(n)
                                lowerTestAct=lowerTestAct+1;
                            end
                        elseif criterionType==2%equal to or higher
                            if cumulativeChAct2(n)<cumulativeChAct4(n)
                                higherTestAct=higherTestAct+1;
                            elseif cumulativeChAct2(n)>=cumulativeChAct4(n)
                                lowerTestAct=lowerTestAct+1;
                            end
                        end
                    end
                end
                rocvals_sglroc3_allCumulative(condInd)=rocMUA(cumulativeChAct4,cumulativeChAct2);%traditional method, can handle negative MUA values
                %rocvals_sglroc3_allCumulative(condInd)=sglroc3(cumulativeChAct4,cumulativeChAct2);%old method, cannot handle negative MUA values
                rocvals_allCumulative(condInd)=higherTestAct/(higherTestAct+lowerTestAct);%new method
                plot(testContrast(condInd),rocvals_sglroc3_allCumulative(condInd),'ro','MarkerSize',3);
                plot(testContrast(condInd),rocvals_allCumulative(condInd),'bo','MarkerSize',3);
            end
            fillCoords=[iqr_roc_upper fliplr(iqr_roc_lower)];
            fillCoords_sglroc3=[iqr_roc_sglroc3_upper fliplr(iqr_roc_sglroc3_lower)];
            xfillCoords=[testContrast fliplr(testContrast)];
            patch(xfillCoords,fillCoords,'b','EdgeColor','none');
            patch(xfillCoords,fillCoords_sglroc3,'r','EdgeColor','none');alpha(0.3)
            xlimVals=get(gca,'XLim');
%             plot([xlimVals(1) xlimVals(2)],[0.5 0.5],'k--');
            if areaInd==2
                xlabel('test contrast (%)','FontAngle','italic');
            end
            if animalInd==1
                ylabel('AUROC or COBAM','FontAngle','italic');
            end
            if calcParams==1
                [slopeNeuro_sglroc3,PNE_sglroc3,diffPNE_sglroc3,minRate_sglroc3,maxRate_sglroc3,chSSE_sglroc3]=weibull_fitting(rocvals_sglroc3_allCumulative,sampleContrast,testContrast,'old',counter,slopeNeuro_sglroc3,[],PNE_sglroc3,minRate_sglroc3,maxRate_sglroc3,diffPNE_sglroc3,0,calculateTangent,0);
                [slopeNeuro,PNE,diffPNE,minRate,maxRate,chSSE]=weibull_fitting(rocvals_allCumulative,sampleContrast,testContrast,'old',counter,slopeNeuro,[],PNE,minRate,maxRate,diffPNE,0,calculateTangent,0);
            end
            xlim([testContrast(1) testContrast(end)]);
            ylim([0 1]);
            axis square
            set(gca, 'box', 'off')
            counter=counter+1;
        end
    end
end
equalCount
notEqualCount
plot2svg%export as an SVG file and open in Illustrator, then save as an EPS/TIF

if calcParams==1
    if calculateTangent==0
        figure;
        subplot(1,2,1);
        plot(slopeNeuro_sglroc3,slopeNeuro,'k.');hold on
        axis square
        ylim([1 6]);
        xlim([1 6]);
        plot([1 6],[1 6],'k--');
        xlabel('AUROC');
        ylabel('COBAM');
        title('Slope');
        subplot(1,2,2);
        range=maxRate-minRate;
        range_sglroc3=maxRate_sglroc3-minRate_sglroc3;
        plot(range_sglroc3,range,'k.');hold on
        axis square
        ylim([0.4 1]);
        plot([0.4 1],[0.4 1],'k--');
        xlabel('AUROC');
        ylabel('COBAM');
        title('Range');
        set(gca,'XTick',[0.4 0.6 0.8 1],'XTickLabel',[0.4 0.6 0.8 1]);
    elseif calculateTangent==1
        figure;
        subplot(1,2,1);
        plot(slopeNeuro_sglroc3,slopeNeuro,'k.');hold on
        axis square
        ylim([0 0.04]);
        xlim([0 0.04]);
        plot([0 0.04],[0 0.04],'k--');
        xlabel('AUROC');
        ylabel('COBAM');
        title('Slope');
        subplot(1,2,2);
        range=maxRate-minRate;
        range_sglroc3=maxRate_sglroc3-minRate_sglroc3;
        plot(range_sglroc3,range,'k.');hold on
        axis square
        ylim([0.4 1]);
        set(gca,'XTick',[0.4 0.6 0.8 1],'XTickLabel',[0.4 0.6 0.8 1]);
        plot([0.4 1],[0.4 1],'k--');
        xlabel('AUROC');
        ylabel('COBAM');
        title('Range');
    end
    [h1,p1,ci1,stats1]=ttest(slopeNeuro,slopeNeuro_sglroc3)
    [h2,p2,ci2,stats2]=ttest(range,range_sglroc3)
end
%use slope param from function:
%t(7) = 4.00, p = .0052
%use slope at 30%:
%t(7) = 5.17, p = .0013

%range:
%t(7) = 2.16, p = .0673