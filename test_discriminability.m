function test_discriminability
%Written by Xing 25/10/13
%Calculates AUROC values for pairs of test contrast, checks for changes
%over the course of training. Uses traditional AUROC approach.
useISI=0;
conciseComparisonConds=1;
excludeSuppressed=0;
cutoff=1;
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
readData=1;
startEndTime=['_1024_to_1536'];
figROCconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
set(figROCconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
figROCcondsCorrIncorr=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
set(figROCcondsCorrIncorr, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
% figROCcondsIncorr=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
% set(figROCcondsIncorr, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
figROCcondsjExample=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
set(figROCcondsjExample, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
if readData==1
    for animalInd=1:length(animals)
        animal=animals{animalInd};
        for areaInd=1:length(areas)
            area=areas{areaInd};
            [sampleContrasts testContrasts]=area_metadata(area);
            if strncmp(area,'v1_2',4)
                sampleContrasts=30;%
                testContrasts=testContrasts(2,:);
            end
            if strncmp(area,'v4_1',4)
                comparisonConds=[5 10;6 9;7 8];%test contrasts of [27 33;28 32;29 31];
                condCol='rgb';
                if conciseComparisonConds==1
                    comparisonConds=[7 8];%test contrasts of [29 31];
                    condCol='b';
                end
            elseif strncmp(area,'v1_1',4)
                comparisonConds=[7 8];%test contrasts of [28 32];
                condCol='g';
            end
            figChconds=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(figChconds, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
            figChcondsCorr=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(figChcondsCorr, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
            figChcondsIncorr=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
            set(figChcondsIncorr, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
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
                corrTrialsAct=cell(length(testContrast),1);
                incorrTrialsAct=cell(length(testContrast),1);
                allEpoch4=cell(length(channels),length(testContrast));
                rocvals_sglroc3=cell(length(channels),1);
                rocvals_sglroc3_corr=cell(length(channels),1);
                rocvals_sglroc3_incorr=cell(length(channels),1);
                sumChActLower=cell(length(sessionNums),size(comparisonConds,1));
                sumChActHigher=cell(length(sessionNums),size(comparisonConds,1));
                sumChActLowerCorr=cell(length(sessionNums),size(comparisonConds,1));
                sumChActHigherCorr=cell(length(sessionNums),size(comparisonConds,1));
                sumChActLowerIncorr=cell(length(sessionNums),size(comparisonConds,1));
                sumChActHigherIncorr=cell(length(sessionNums),size(comparisonConds,1));
                rhoChCorr=[];
                pSpearmanChCorr=[];
                dfsChCorr=[];
                rhoChIncorr=[];
                pSpearmanChIncorr=[];
                dfsChIncorr=[];
                for chInd=1:length(channels)
                    for i=1:length(sessionNums)%combine act across trials and across all sessions for each channel
                        for condInd=1:size(comparisonConds,1)
                            matPath=['F:\PL\sample_test_activity\',animal,'_',area,'\ch',num2str(channels(chInd)),'_',num2str(sessionNums(i)),'_example_sample_test_act.mat'];
                            loadText=['load ',matPath];
                            eval(loadText);
                            corr_incorr_matName=[area,'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'.mat'];
                            corr_incorr_matPath=fullfile('F:','PL','TAS',animal,corr_incorr_matName);
                            loadText=['load ',corr_incorr_matPath,' respArray'];
                            eval(loadText);
                            corrInd=find(respArray{comparisonConds(condInd,2)}==1);
                            incorrInd=find(respArray{comparisonConds(condInd,2)}==-1);
                            corrTrialsAct{comparisonConds(condInd,2)}=epoch4{comparisonConds(condInd,2)}(corrInd);
                            incorrTrialsAct{comparisonConds(condInd,2)}=epoch4{comparisonConds(condInd,2)}(incorrInd);
                            corrInd=find(respArray{comparisonConds(condInd,1)}==1);
                            incorrInd=find(respArray{comparisonConds(condInd,1)}==-1);
                            corrTrialsAct{comparisonConds(condInd,1)}=epoch4{comparisonConds(condInd,1)}(corrInd);
                            incorrTrialsAct{comparisonConds(condInd,1)}=epoch4{comparisonConds(condInd,1)}(incorrInd);
                            matExists=0;
                            if exist(matPath,'file')
                                matExists=1;
                            end
                            includeRows=includeSessionsAll(find(includeSessionsAll(:,1)==channels(chInd)),2);%include this session in analysis
                            includeRow=find(includeRows==sessionNums(i));
                            eval(loadText);
                            %calculate ROC values for each channel and session:
                            rocvals_sglroc3{chInd}(i,condInd)=sglroc3(epoch4{comparisonConds(condInd,2)},epoch4{comparisonConds(condInd,1)});
                            rocvals_sglroc3_corr{chInd}(i,condInd)=sglroc3(corrTrialsAct{comparisonConds(condInd,2)},corrTrialsAct{comparisonConds(condInd,1)});
                            rocvals_sglroc3_incorr{chInd}(i,condInd)=sglroc3(incorrTrialsAct{comparisonConds(condInd,2)},incorrTrialsAct{comparisonConds(condInd,1)});
                            %combine activity across channels for each trial:
                            if isempty(sumChActLower{i,condInd})
                                sumChActLower{i,condInd}=epoch4{comparisonConds(condInd,1)};
                                sumChActLowerCorr{i,condInd}=corrTrialsAct{comparisonConds(condInd,1)};
                                sumChActLowerIncorr{i,condInd}=incorrTrialsAct{comparisonConds(condInd,1)};
                            else
                                sumChActLower{i,condInd}=sumChActLower{i,condInd}+epoch4{comparisonConds(condInd,1)};
                                sumChActLowerCorr{i,condInd}=sumChActLowerCorr{i,condInd}+corrTrialsAct{comparisonConds(condInd,1)};
                                sumChActLowerIncorr{i,condInd}=sumChActLowerIncorr{i,condInd}+incorrTrialsAct{comparisonConds(condInd,1)};
                            end
                            if isempty(sumChActHigher{i,condInd})
                                sumChActHigher{i,condInd}=epoch4{comparisonConds(condInd,2)};
                                sumChActHigherCorr{i,condInd}=corrTrialsAct{comparisonConds(condInd,2)};
                                sumChActHigherIncorr{i,condInd}=incorrTrialsAct{comparisonConds(condInd,2)};
                            else
                                sumChActHigher{i,condInd}=sumChActHigher{i,condInd}+epoch4{comparisonConds(condInd,2)};
                                sumChActHigherCorr{i,condInd}=sumChActHigherCorr{i,condInd}+corrTrialsAct{comparisonConds(condInd,2)};
                                sumChActHigherIncorr{i,condInd}=sumChActHigherIncorr{i,condInd}+incorrTrialsAct{comparisonConds(condInd,2)};
                            end
                        end
                    end
                    figure(figChconds);
                    subplot(ceil(length(channels)/5),5,chInd);
                    for condInd=1:size(comparisonConds,1)
                        plot(1:length(sessionNums),rocvals_sglroc3{chInd}(:,condInd),'Color',condCol(condInd),'Marker','o');hold on
                    end
                    title(num2str(channels(chInd)));
                    if chInd==1
                        xlabel('session number','FontAngle','italic');
                        ylabel('AUROC value (all trials)','FontAngle','italic');
                    end
                    figure(figChcondsCorr);
                    subplot(ceil(length(channels)/5),5,chInd);
                    for condInd=1:size(comparisonConds,1)
                        plot(1:length(sessionNums),rocvals_sglroc3_corr{chInd}(:,condInd),'Color',condCol(condInd),'Marker','o');hold on
                    end
                    title(num2str(channels(chInd)));
                    if chInd==1
                        xlabel('session number','FontAngle','italic');
                        ylabel('AUROC value','FontAngle','italic');
                    end
                    figure(figChcondsIncorr);
                    subplot(ceil(length(channels)/5),5,chInd);
                    for condInd=1:size(comparisonConds,1)
                        plot(1:length(sessionNums),rocvals_sglroc3_incorr{chInd}(:,condInd),'Color',condCol(condInd),'Marker','o');hold on
                    end
                    title(num2str(channels(chInd)));
                    if chInd==1
                        xlabel('session number','FontAngle','italic');
                        ylabel('AUROC value (incorrect trials)','FontAngle','italic');
                    end
                    
                    if channels(chInd)==53&&strcmp(animal,'jack')||channels(chInd)==55&&strcmp(animal,'jack')
                        figure(figROCcondsjExample)
                        if channels(chInd)==53
                            subplot(1,2,1);
                        else
                            subplot(1,2,2);
                        end
                        for condInd=1:size(comparisonConds,1)
                            plot(1:length(sessionNums),rocvals_sglroc3_corr{chInd}(:,condInd),'Color',condCol(condInd),'Marker','o','MarkerFaceColor',condCol(condInd),'LineStyle','none');hold on
                            plot(1:length(sessionNums),rocvals_sglroc3_incorr{chInd}(:,condInd),'Color',condCol(condInd),'Marker','o','LineStyle','none');hold on
                            hold on
                        end
                        title(num2str(channels(chInd)));
                        if channels(chInd)==53
                            ylabel('AUROC value','FontAngle','italic');
                        end
                        xlabel('session number','FontAngle','italic');
                        axis square
                        set(gca, 'box', 'off')
                        xlimVals=get(gca,'XLim');
                        plot([xlimVals(1) xlimVals(2)],[0.5 0.5],'k--');                        
                    end
                    
                    %calculate stats:
                    c=[(1:length(sessionNums))' (rocvals_sglroc3_corr{chInd}(:,condInd))];%session in column 1, data in column 2
                    [a b]=corr(c,'type','Spearman');
                    rhoChCorr(chInd)=a(2);
                    pSpearmanChCorr(chInd)=b(2);
                    dfsChCorr(chInd)=size(c,1)-2;
                    c=[(1:length(sessionNums))' (rocvals_sglroc3_incorr{chInd}(:,condInd))];%session in column 1, data in column 2
                    [a b]=corr(c,'type','Spearman');
                    rhoChIncorr(chInd)=a(2);
                    pSpearmanChIncorr(chInd)=b(2);
                    dfsChIncorr(chInd)=size(c,1)-2;
                end
                figure(figROCconds);
                subplot(2,2,animalInd+2*(areaInd-1));
                sum_rocvals_sglroc3=[];
                for condInd=1:size(comparisonConds,1)
                    for i=1:length(sessionNums)%combine act across trials and across all sessions for each channel
                        sum_rocvals_sglroc3(i,condInd)=sglroc3(sumChActHigher{i,condInd},sumChActLower{i,condInd});%activity summed across channels on each trial
                    end
                    plot(1:size(sum_rocvals_sglroc3,1),sum_rocvals_sglroc3(:,condInd),'Color',condCol(condInd),'Marker','o');hold on
                end
                if areaInd==2
                    xlabel('session number','FontAngle','italic');
                end
                if animalInd==1
                    ylabel('AUROC value (all trials)','FontAngle','italic');
                end
                set(gca, 'box', 'off')
                figure(figROCcondsCorrIncorr);
                subplot(2,2,animalInd+2*(areaInd-1));
                sum_rocvals_sglroc3_corr=[];
                sum_rocvals_sglroc3_incorr=[];
                for condInd=1:size(comparisonConds,1)
                    for i=1:length(sessionNums)%combine act across trials and across all sessions for each channel
                        sum_rocvals_sglroc3_corr(i,condInd)=sglroc3(sumChActHigherCorr{i,condInd},sumChActLowerCorr{i,condInd});%activity summed across channels on each trial
                        sum_rocvals_sglroc3_incorr(i,condInd)=sglroc3(sumChActHigherIncorr{i,condInd},sumChActLowerIncorr{i,condInd});%activity summed across channels on each trial
                    end
                    plot(1:size(sum_rocvals_sglroc3_corr,1),sum_rocvals_sglroc3_corr(:,condInd),'Color',condCol(condInd),'Marker','o','MarkerFaceColor',condCol(condInd),'LineStyle','none');hold on
                    plot(1:size(sum_rocvals_sglroc3_incorr,1),sum_rocvals_sglroc3_incorr(:,condInd),'Color',condCol(condInd),'Marker','o','MarkerFaceColor','none','LineStyle','none');hold on
                end
                if areaInd==2
                    xlabel('session number','FontAngle','italic');
                end
                if animalInd==1
                    ylabel('AUROC value','FontAngle','italic');
                end
                set(gca, 'box', 'off')
                xlimVals=get(gca,'XLim');
                plot([xlimVals(1) xlimVals(2)],[0.5 0.5],'k--');
                
                %calculate stats:
                c=[(1:size(sum_rocvals_sglroc3_corr,1))' (sum_rocvals_sglroc3_corr(:,condInd))];%session in column 1, data in column 2
                [a b]=corr(c,'type','Spearman');
                rhoCorr(animalInd+2*(areaInd-1))=a(2);
                pSpearmanCorr(animalInd+2*(areaInd-1))=b(2);
                dfsCorr(animalInd+2*(areaInd-1))=size(c,1)-2;
                c=[(1:size(sum_rocvals_sglroc3_incorr,1))' (sum_rocvals_sglroc3_incorr(:,condInd))];%session in column 1, data in column 2
                [a b]=corr(c,'type','Spearman');
                rhoIncorr(animalInd+2*(areaInd-1))=a(2);
                pSpearmanIncorr(animalInd+2*(areaInd-1))=b(2);
                dfsIncorr(animalInd+2*(areaInd-1))=size(c,1)-2;
%                 figure(figROCcondsIncorr);
%                 subplot(2,2,animalInd+2*(areaInd-1));
%                 sum_rocvals_sglroc3_incorr=[];
%                 for condInd=1:size(comparisonConds,1)
%                     for i=1:length(sessionNums)%combine act across trials and across all sessions for each channel
%                         sum_rocvals_sglroc3_incorr(i,condInd)=sglroc3(sumChActHigherIncorr{i,condInd},sumChActLowerIncorr{i,condInd});%activity summed across channels on each trial
%                     end
%                     plot(1:size(sum_rocvals_sglroc3,1),sum_rocvals_sglroc3_incorr(:,condInd),'Color',condCol(condInd),'Marker','o');hold on
%                 end
%                 if areaInd==2
%                     xlabel('session number','FontAngle','italic');
%                 end
%                 if animalInd==1
%                     ylabel('AUROC value (incorrect trials)','FontAngle','italic');
%                 end
%                 set(gca, 'box', 'off')
            end
        end
    end
end
     