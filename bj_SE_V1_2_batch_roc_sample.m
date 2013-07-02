function bj_SE_V1_2_batch_roc_sample
%Modified from blanco_SE_V1_2_batch_roc_sample on 26/03/13 to analyse
%activity during sample presentation and calculate ROC values between
%firing to 20 and 30% sample, and between 30 and 40% sample.

animalTexts=[{'Monkey 1'} {'Monkey 2'}];
animals=[{'blanco'} {'jack'}];
areaText='V1 roving data';
areas={'v1_2_1' 'v1_2_2' 'v1_2_3'};
readSampleAct=0;
plotIndividualChs=0;
markerCols='kcm';
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        sigList=[];
        sigCondList=[];
        if nargin<3 || isempty(channels)
            channels = main_channels(animal,area);
        end
        if nargin<4 || isempty(sessionNums)
            sessionNums = main_raw_sessions_final(animal,area,[],0);
        end
        [sampleContrasts allTestContrasts]=area_metadata(area);
        epoch=2;%calculate ROC values based on response to sample stimuli
        minusSpontan=0;
        if minusSpontan==1
            subfolder=['roving_sample_ROC',num2str(epoch),'_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
        elseif minusSpontan==0
            subfolder=['roving_sample_ROC',num2str(epoch),'_wspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
        end
        if readSampleAct==1
            for i=1:length(channels)
                for j=1:length(sessionNums)
                    for sampleContrastInd=1:length(sampleContrasts)
                        sampleContrast=sampleContrasts(sampleContrastInd);
                        testContrasts=allTestContrasts(sampleContrastInd,:);
                        try
                            bj_SE_roving_sample_ROC(animal,area,channels(i),sessionNums(j),sampleContrast,testContrasts,epoch);%note that epochTimes(ind,2) is currently unused
                        catch ME
                            disp(ME)
                            load F:\PL\roving_sample_ROC\missingSessions.mat missingSessions
                            missingSessions=[missingSessions;{animal} {channels(i)} {sessionNums(j)} {ME}];
                            save F:\PL\roving_sample_ROC\missingSessions.mat missingSessions
                        end
                    end
                end
            end
        end
        if plotIndividualChs==1
            figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
            for i=1:length(channels)
                subplot(ceil(length(channels)/5),5,i);
                rocvals=[];
                allrocvals=[];
                for j=1:length(sessionNums)
                    for sampleContrastInd=1:length(sampleContrasts)
                        sampleContrast=sampleContrasts(sampleContrastInd);
                        matSampleName=[num2str(channels(i)),'_',num2str(sessionNums(j)),'_',num2str(sampleContrast),'_',area,'_sample_vals'];
                        matSampleFolder=fullfile('F:','PL','roving_sample_ROC',animal,subfolder);
                        if ~exist(matSampleFolder,'dir')
                            mkdir(matSampleFolder)
                        end
                        matSamplePath=fullfile(matSampleFolder,matSampleName);
                        loadText=['load ',matSamplePath,' sample_act'];
                        eval(loadText);
                        allSamples(j,sampleContrastInd)={sample_act};%compile across 3 sample contrast conditions
                    end
                    %roc analysis:
                    compareSamples=[2 1;3 2;3 1];
                    for comparisonInd=1:length(sampleContrasts)
                        [roc]=sglroc3(allSamples{j,compareSamples(comparisonInd,1)}',allSamples{j,compareSamples(comparisonInd,2)}');
                        rocvals(1,comparisonInd)=roc;%20 vs 30; 30 vs 40; 20 vs 40
                    end
                    allrocvals(j,:)=rocvals;
                end
                sigChFlag=0;
                for comparisonInd=1:length(sampleContrasts)
                    a=[(1:size(allrocvals,1))' allrocvals(:,comparisonInd)];
                    [coefficients1 p1]=corrcoef(a);
                    coefficients(comparisonInd)=coefficients1(2);
                    ps(comparisonInd)=p1(2);
                    if p1(2)<0.05/3
                        markerCol2='k';
                        sigCondList=[sigCondList;channels(i) sampleContrasts(comparisonInd) p1(2) coefficients1(2)];
                        sigChFlag=1;
                    else
                        markerCol2='none';
                    end
                    plot(1:length(sessionNums),allrocvals(:,comparisonInd),'Color',markerCols(comparisonInd),'Marker','o','MarkerFaceColor',markerCol2,'LineStyle','none');hold on
                end
                if sigChFlag==1
                    sigList=[sigList;channels(i)];                    
                end
                if i==length(channels)
                    xlabel('session number');
                elseif i==1
                    ylabel('ROC');
                end
                minY=min(min(allrocvals))-0.05;
                maxY=max(max(allrocvals));
                ylim([minY maxY]);
                ylims=get(gca,'yLim');
                ptext=sprintf('R= %.3f %.3f %.3f p= %.3f %.3f %.3f',coefficients(1),coefficients(2),coefficients(3),ps(1),ps(2),ps(3));
                text('Position',[0 0.05*(maxY-minY)+minY],'FontSize',9,'String',ptext);
                title(channels(i));        %plot ROC data points for all channels (not averaged) against sessions number
                sampleROCs(i,:)=[channels(i) {allSamples} {allrocvals} {coefficients} {ps}];
            end
            saveImageName=[area,'_',animal,'_all_ch_sample_ROC'];
            saveImageFolder=fullfile('F:','PL','roving_sample_ROC',animal);
            saveImagePath=fullfile(saveImageFolder,saveImageName);
            printtext=sprintf('print -dpng -r600 %s.png',saveImagePath);
            set(gcf,'PaperPositionMode','auto')
            eval(printtext);
            
            matSampleName=[animal,'_',area,'_all_ch_sample_ROC'];
            matSampleFolder=fullfile('F:','PL','roving_sample_ROC',animal,subfolder);
            if ~exist(matSampleFolder,'dir')
                mkdir(matSampleFolder)
            end
            matSamplePath=fullfile(matSampleFolder,matSampleName);
            saveText=['save ',matSamplePath,' sampleROCs'];
            eval(saveText);
        end
    end
end

areaTexts=[{'no flankers'} {'flankers'}];
figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.1, 0.1, 0.4, 0.5]);
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:2%exclude the final V1_3 phase in which flankers are removed
        area=areas{areaInd};
        channels = main_channels(animal,area);
        sessionNums = main_raw_sessions_final(animal,area,[],0);
        subplot(2,2,animalInd+2*(areaInd-1));
        if minusSpontan==1
            subfolder=['roving_sample_ROC',num2str(epoch),'_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
        elseif minusSpontan==0
            subfolder=['roving_sample_ROC',num2str(epoch),'_wspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
        end
        matSampleName=[animal,'_',area,'_all_ch_sample_ROC'];
        matSampleFolder=fullfile('F:','PL','roving_sample_ROC',animal,subfolder);
        if ~exist(matSampleFolder,'dir')
            mkdir(matSampleFolder)
        end
        matSamplePath=fullfile(matSampleFolder,matSampleName);
        loadText=['load ',matSamplePath,' sampleROCs'];
        eval(loadText);
        for comparisonInd=1:3
            meanSampleROCs=zeros(length(sessionNums),1);
            for rowInd=1:length(channels)
                meanSampleROCs=meanSampleROCs+sampleROCs{rowInd,3}(:,comparisonInd);
            end
            meanSampleROCs=meanSampleROCs/length(channels);
            d=[(1:length(meanSampleROCs))' meanSampleROCs];
            [coefficients1 p1]=corrcoef(d);
            coefficients(comparisonInd)=coefficients1(2);
            ps(comparisonInd)=p1(2);
            if p1(2)<0.05/3
                markerCol2=markerCols(comparisonInd);
            else
                markerCol2='none';
            end
            plot(1:length(sessionNums),meanSampleROCs,'Color',markerCols(comparisonInd),'Marker','o','MarkerFaceColor',markerCol2,'LineStyle','none');
            hold on
            allMeanSampleROCs{comparisonInd}=meanSampleROCs;
        end
        xlim([1 length(meanSampleROCs)]);
        ylims=get(gca,'YLim');
        set(gca,'YTick',[ylims(1) ylims(2)]);
        set(gca,'YTickLabel',[ylims(1) ylims(2)]);
        if animalInd==1
            ylabel('ROC value');
            if areaInd==2
                text('Position',[-10 0.5],'FontSize',9,'String',areaTexts{areaInd});
            end
        end
        if areaInd==1
            title(animalTexts{animalInd});
        else
            xlabel('session');
        end
        matMeanSampleName=[animal,'_',area,'_mean_ch_sample_ROC'];
        matMeanSampleFolder=fullfile('F:','PL','roving_sample_ROC',animal,subfolder);
        if ~exist(matMeanSampleFolder,'dir')
            mkdir(matMeanSampleFolder)
        end
        matSamplePath=fullfile(matMeanSampleFolder,matMeanSampleName);
        saveText=['save ',matSamplePath,' allMeanSampleROCs coefficients ps'];
        eval(saveText);
    end
end
subplot(2,2,1)
set(gca,'YTick',[0.45 0.6]);
ylim([0.45 0.6]);
subplot(2,2,3)
set(gca,'YTick',[0.45 0.57],'YTickLabel',[0.45 0.57]);
ylim([0.45 0.57]);
saveMeanImageName=[area,'_mean_ch_sample_ROC'];
saveMeanImageFolder=fullfile('F:','PL','roving_sample_ROC');
saveMeanImagePath=fullfile(saveMeanImageFolder,saveMeanImageName);
printtext=sprintf('print -dpng -r600 %s.png',saveMeanImagePath);
% axis square
set(gcf,'PaperPositionMode','auto')
eval(printtext);
