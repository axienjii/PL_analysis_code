function bj_roc_diff_figure_vs_corr(exampleFig,cutoff,animals,useISI,areas,excludeSuppressed,normalize,paperSessions,readMehdiSNR)
%Written by Xing 28/05/14
%Modified from bj_roc_diff_figure, run current function after reading data 
%from this previous function.
%Generate figure for methods paper, comparing AUROC and COBAM methods.
%Read activity levels for all trials (separately for each session and test
%contrast) across channels, calculate population AUROC and DICAF values.
%Flip DICAF and AUROC values about 0.5 as needed.
%Calculate within-trial sample-test activity correlation, R, for each session,
%plot difference index (DICAF-AUROC)/AUROC against R.

%Set useISI to 1: based on pre-test vs test, not on sample vs test.
%Set useISI to 0: sample vs test.

%Set exampleFig to 0 to plot ROC curves for new and old methods, set to 1
%to only plot example figures of distributions of stimulus-evoked activity
%and condition-dependent ROC curves.
%Set excludeSuppressed to 1 to exclude channels with stimulus-evoked
%suppression, i.e. blanco 13, 24, 42 and jack 49.
plotRsSeparateAnimals=0;
criterionType=0;
calculateTangent=1;
equalCount=0;
notEqualCount=0;
calcParams=1;
if useISI==1
    analysisType='ROC';
else
    analysisType='ROC_zero_one';
end
counter=1;%8 values- spiking B V4, B V1, J V4, J V1, MUA B V4, B V1, J V4, J V1
animalCol=[{[0.4 0.4 0.4]} {'k'}];
areaMarker=[{'none'} {'k'}];
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
readData=1;
figSessRs=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %AUROC & DICAF vals against within-trial Rs
set(figSessRs, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
figDIsRs=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %difference indices against within-trial Rs
set(figDIsRs, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
figDIsRsAll=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.3, 0.2]);
set(figDIsRsAll, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
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
            mixedDirections=0;
            if readMehdiSNR==1%read in list of included channels
                if cutoff~=1
                    matname=['good_SNR_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10),'.mat'];
                    pathname=fullfile(rootFolder,'PL','SNR',animal,matname);
                else
                    matname=['good_SNR_',area,'_',num2str(sampleContrast),'.mat'];
                    pathname=fullfile(rootFolder,'PL','SNR',animal,'cutoff_SNR_1',matname);
                end
                loadText=['load ',pathname,' includeSessionsAll'];
                eval(loadText);
            end
            rocvals_sglroc3=zeros(length(testContrast),length(sessionNums));
            rocvals=zeros(length(testContrast),length(sessionNums));
            DI=zeros(length(testContrast),length(sessionNums));
            withinTrialR=zeros(length(testContrast),length(sessionNums));
            if useISI==0
                allEpoch2=cell(length(channels),length(testContrast));%array to store activity for each ch and cond across sessions and trials
            elseif useISI==1
                allEpoch3=cell(length(channels),length(testContrast));
            end
            epoch2Sessions=cell(length(testContrast),length(sessionNums));%sum act across chs for each test contrast and session. Rows: test contrasts; columns: sessions
            epoch4Sessions=cell(length(testContrast),length(sessionNums));
            for i=1:length(sessionNums)%combine act across channels for each trial
                for chInd=1:length(channels)
                    matPath=['F:\PL\sample_test_activity\',animal,'_',area,'\ch',num2str(channels(chInd)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat'];
                    loadText=['load ',matPath];
                    matExists=0;
                    if exist(matPath,'file')
                        matExists=1;
                    end
                    if readMehdiSNR
                        includeRows=includeSessionsAll(find(includeSessionsAll(:,1)==channels(chInd)),2);%include this session in analysis
                        includeRow=find(includeRows==sessionNums(i));
                        SNRgood=~isempty(includeRow);
                    else
                        SNRgood=0;
                    end
                    if matExists==1&&SNRgood==1||matExists==1&&readMehdiSNR==0%&&includeCh
                        eval(loadText);
                        for condInd=1:length(testContrast)
                            if useISI==0
                                if isempty(epoch2Sessions{condInd,i})
                                    epoch2Sessions{condInd,i}=epoch2{condInd};
                                else
                                    epoch2Sessions{condInd,i}=[epoch2Sessions{condInd,i}+epoch2{condInd}];
                                end
                            elseif useISI==1
                                if isempty(epoch2Sessions{condInd,i})
                                    epoch3Sessions{condInd,i}=epoch3{condInd};
                                else
                                    epoch3Sessions{condInd,i}=[epoch3Sessions{condInd,i}+epoch3{condInd}];
                                end
                            end
                            if isempty(epoch4Sessions{condInd,i})
                                epoch4Sessions{condInd,i}=epoch4{condInd};
                            else
                                epoch4Sessions{condInd,i}=[epoch4Sessions{condInd,i}+epoch4{condInd}];
                            end
                            if normalize==1
                                maxAll=[];%find the highest firing rate across all conditions and trials, across both the sample and test presentation periods
                                maxval=max(maxAll)/100;
                            else
                                maxval=1;
                            end
                        end
                        if chInd==1
                            numSessTrials(i,condInd)=size(epoch2{condInd},2);
                        end
                        %                             if chInd==15
                        %                                 numSessTrials49(i,condInd)=size(epoch2{condInd},2);
                        %                             end
                    end
                end
                %for each session and test contrast condition, calculate R, AUROC, DICAF, and DI
                for condInd=1:length(testContrast)
                    higherTestAct=0;
                    lowerTestAct=0;
                    for n=1:length(epoch2Sessions{condInd,i})%tally across trials
                        if useISI==0
                            if criterionType==0%exclude trials where responses to sample and test are equal
                                if epoch2Sessions{condInd,i}(n)<epoch4Sessions{condInd,i}(n)
                                    higherTestAct=higherTestAct+1;
                                elseif epoch2Sessions{condInd,i}(n)>epoch4Sessions{condInd,i}(n)
                                    lowerTestAct=lowerTestAct+1;
                                end
                            elseif criterionType==1%for trials where responses to sample and test are equal, assign to 'higher test' category
                                if epoch2Sessions{condInd,i}(n)<=epoch4Sessions{condInd,i}(n)
                                    higherTestAct=higherTestAct+1;
                                elseif epoch2Sessions{condInd,i}(n)>epoch4Sessions{condInd,i}(n)
                                    lowerTestAct=lowerTestAct+1;
                                end
                            elseif criterionType==2%for trials where responses to sample and test are equal, assign to 'lower test' category
                                if epoch2Sessions{condInd,i}(n)<epoch4Sessions{condInd,i}(n)
                                    higherTestAct=higherTestAct+1;
                                elseif epoch2Sessions{condInd,i}(n)>=epoch4Sessions{condInd,i}(n)
                                    lowerTestAct=lowerTestAct+1;
                                end
                            end
                            if epoch2Sessions{condInd,i}(n)==epoch4Sessions{condInd,i}(n)
                                equalCount=equalCount+1;
                            else
                                notEqualCount=notEqualCount+1;
                            end
                        elseif useISI==1
                            if epoch3Sessions{condInd,i}(n)<epoch4Sessions{condInd,i}(n)
                                higherTestAct=higherTestAct+1;
                            elseif epoch3Sessions{condInd,i}(n)>epoch4Sessions{condInd,i}(n)
                                lowerTestAct=lowerTestAct+1;
                            end
                        end
                    end
                    if isempty(epoch4Sessions{condInd,i})||isempty(epoch2Sessions{condInd,i})
                        pauseHere=1;
                    else
                        rocvals_sglroc3(condInd,i)=sglroc3(epoch4Sessions{condInd,i},epoch2Sessions{condInd,i});%old method
                        rocvals(condInd,i)=higherTestAct/(higherTestAct+lowerTestAct);%new method
                        [R p]=corrcoef([epoch4Sessions{condInd,i}' epoch2Sessions{condInd,i}']);
                        withinTrialR(condInd,i)=R(2);
                        withinTrialp(condInd,i)=p(2);
                        figure(figSessRs)
                        subplot(2,2,animalInd+2*(areaInd-1));
                        if rocvals_sglroc3(condInd,i)>0.5&&rocvals(condInd,i)>0.5
                            DI(condInd,i)=(rocvals(condInd,i)-rocvals_sglroc3(condInd,i))/rocvals_sglroc3(condInd,i);%difference index=(DICAF-AUROC)/AUROC
                            plot(withinTrialR(condInd,i),rocvals_sglroc3(condInd,i),'ro');hold on
                            plot(withinTrialR(condInd,i),rocvals(condInd,i),'bo');
                        elseif rocvals_sglroc3(condInd,i)<0.5&&rocvals(condInd,i)<0.5
                            DI(condInd,i)=(rocvals_sglroc3(condInd,i)-rocvals(condInd,i))/(1-rocvals_sglroc3(condInd,i));%difference index=(AUROC-DICAF)/flipped AUROC
                            plot(withinTrialR(condInd,i),1-rocvals_sglroc3(condInd,i),'ro');hold on
                            plot(withinTrialR(condInd,i),1-rocvals(condInd,i),'bo');
                        else
                            mixedDirections=mixedDirections+1;%AUROC and DICAF values are on opposite sides of 0.5, weird!
                            if rocvals_sglroc3(condInd,i)>0.5&&rocvals(condInd,i)<0.5
                                DI(condInd,i)=((1-rocvals(condInd,i))-rocvals_sglroc3(condInd,i))/rocvals_sglroc3(condInd,i);%difference index=(flipped DICAF-AUROC)/AUROC
                                plot(withinTrialR(condInd,i),rocvals_sglroc3(condInd,i),'ro');hold on
                                plot(withinTrialR(condInd,i),1-rocvals(condInd,i),'bo');
                            elseif rocvals_sglroc3(condInd,i)<0.5&&rocvals(condInd,i)>0.5
                                DI(condInd,i)=(rocvals(condInd,i)-(1-rocvals_sglroc3(condInd,i)))/(1-rocvals_sglroc3(condInd,i));%difference index=(DICAF-flipped AUROC)/flipped AUROC
                                plot(withinTrialR(condInd,i),1-rocvals_sglroc3(condInd,i),'ro');hold on
                                plot(withinTrialR(condInd,i),rocvals(condInd,i),'bo');
                            end
                        end
                        %regardless of the relationships between DICAF,
                        %AUROC, and 0.5, as long as DICAF is further from
                        %0.5 than is AUROC, DI will be positive. If AUROC
                        %is further from 0.5 than DICAF, then DI will be
                        %negative.
                    end
                end
            end
            if plotRsSeparateAnimals==1
                for i=1:length(sessionNums)%combine act across channels for each trial
                    for condInd=1:length(testContrast)
                        figure(figDIsRs)
                        subplot(2,2,animalInd+2*(areaInd-1));
                        plot(withinTrialR(condInd,i),DI(condInd,i),'ko');%plots of DI against within-trial R, based on activity summed across all channels
                        hold on
                    end
                end
            end
            allsessROCvals_sglroc3{areaInd,animalInd}=rocvals_sglroc3;
            allsessROCvals{areaInd,animalInd}=rocvals;
            allRvals{areaInd,animalInd}=withinTrialR;
            allpvals{areaInd,animalInd}=withinTrialp;
            allDIvals{areaInd,animalInd}=DI;
            allSubjectsEpoch4{areaInd,animalInd}=epoch4Sessions;
            allSubjectsEpoch2{areaInd,animalInd}=epoch2Sessions;
            allChannels{areaInd,animalInd}=channels;
            allMixedDirections{areaInd,animalInd}=mixedDirections;
            figROCconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(figROCconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
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
        if paperSessions
            subFolder=[subFolder,'_paperSessions'];
        end
        if useISI==0
            imagename=['population_discrim_vs_withintrialR_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
        elseif useISI==1
            imagename=['population_discrim_vs_withintrialR_old_new_useISI_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
        end
        folderpathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder);
        if ~exist(folderpathname,'dir')
            mkdir(folderpathname);
        end
        pathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder,imagename);
        printtext=sprintf('print -dpng %s.png',pathname);
        set(gcf,'PaperPositionMode','auto')
        eval(printtext);
        %plot DIs against Rs for all subjects & areas
        figure(figDIsRsAll)
        subplot(1,2,areaInd);%data from all subjects on same plot, but separate plots for each area
        for i=1:length(sessionNums)%combine act across channels for each trial
            for condInd=1:length(testContrast)
                plot(withinTrialR,DI,'Marker','o','MarkerFaceColor','none','Color',animalCol{animalInd},'LineStyle','none','MarkerSize',2);%plots of DI against within-trial R, based on activity summed across all channels  %areaMarker{animalInd}
                hold on
            end
        end
        ylimvals=get(gca,'YLim');
        xlimvals=get(gca,'XLim');
        plot([0 0],ylimvals,'k--');
        plot(xlimvals,[0 0],'k--');
        [DIvsRsGrandR DIvsRsGrandp]=corrcoef(withinTrialR(:),DI(:));%calculate correlation between difference indices and within-trial correlation values
        allDIvsRsGrandR{areaInd,animalInd}=DIvsRsGrandR(2);%store across monkeys and areas
        allDIvsRsGrandp{areaInd,animalInd}=DIvsRsGrandp(2);
        ylim([-0.2 0.3]);
        set(gca, 'box', 'off')
    end
end
%allDIvsRsGrandR:
%       B           J
% V4   [0.4138]    [0.2881]
% V1   [0.4197]    [0.2168]

%allDIvsRsGrandp:
%       B           J
% V4   [1.1985e-017]    [4.0869e-008]
% V1   [1.4174e-011]    [1.2562e-004]

%df:
%       B      J
% V4   392    350
% V1   238    308
    
%vertically concatenate all DI and R vals across monkeys & areas, calculate correlation between DIs & within-trial Rs for data pooled across monkeys & areas:
vertDIvals=[allDIvals{1,1}(:);allDIvals{1,2}(:);allDIvals{2,1}(:);allDIvals{2,2}(:)];
vertRvals=[allRvals{1,1}(:);allRvals{1,2}(:);allRvals{2,1}(:);allRvals{2,2}(:)];
[allR allp]=corrcoef([vertDIvals vertRvals]);
grandDIsRsR=allR(2);%0.4112
grandDIsRsp=allp(2);%1.0473e-053

for animalInd=1:length(animals)
    for areaInd=1:length(areas)%calculate mean and SD of within-trial Rs for each animal and area
        meanAllRvals{areaInd,animalInd}=mean(allRvals{areaInd,animalInd}(:));
        SDAllRvals{areaInd,animalInd}=std(allRvals{areaInd,animalInd}(:));
    end
    [h p ci stats]=ttest2(allRvals{1,animalInd}(:),allRvals{2,animalInd}(:))%unpaired t-test, check whether within-trial Rs are significantly different between areas for each animal
end

%meanAllRvals: 
%       B           J
% V4    [0.2260]    [0.1527]
% V1    [0.4343]    [0.2367]

%SDAllRvals: 
%       B           J
% V4    [0.1718]    [0.1860]
% V1    [0.1623]    [0.1515]

% monkey 1
% p: 0
% ci:
%    -0.2354
%    -0.1811
%    stats:
%     tstat: -15.0634
%        df: 628
%        sd: 0.1682

%monkey 2
% p: 5.7120e-010
% ci:
%    -0.1101
%    -0.0577
%    stats:
%     tstat: -6.2925
%        df: 656
%        sd: 0.1707

if excludeSuppressed==0
    saveText=['save ',rootFolder,'\PL\',analysisType,'\population_AUROC_DICAF_Rvals allsessROCvals_sglroc3 allsessROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels allDIvsRsGrandR allDIvsRsGrandp allRvals allpvals grandDIsRsR grandDIsRsp allMixedDirections meanAllRvals SDAllRvals'];
    eval(saveText)
elseif excludeSuppressed==1
    saveText=['save ',rootFolder,'\PL\',analysisType,'\population_AUROC_DICAF_Rvals_excludeSuppressed allsessROCvals_sglroc3 allsessROCvals allSubjectsEpoch2 allSubjectsEpoch4 allChannels allDIvsRsGrandR allDIvsRsGrandp allRvals allpvals grandDIsRsR grandDIsRsp allMixedDirections meanAllRvals SDAllRvals'];
    eval(saveText)
end