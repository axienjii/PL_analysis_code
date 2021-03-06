function bj_individual_ch_new_old_roc_diff(exampleFig,cutoff,animals,useISI,areas,excludeSuppressed,normalize)
%Written by Xing 24/07/13- the analyses carried out by this function were
%previously done in bj_cumulative_roc_diff, until that function got
%modified to analyse population activity rather than single channel
%activity. Thus current function written to take over that role.

%Set useISI to 1: based on pre-test vs test, not on sample vs test.
%Set useISI to 0: sample vs test.
%Calculate ROC values for individual channels.
%To compare cumulatively-calculated ROC values obtained using 2 methods:
%1. 'AUROC' sglroc3 (blue, old) and 2. 'PROBMAT' the mean trial-wise higher/lower activity (red, new).
%Also plots distributions of sample- and test-evoked activity, and ROC
%curves, for each of the two methods.
%Set exampleFig to 0 to plot ROC curves for new and old methods, set to 1
%to only plot example figures of distributions of stimulus-evoked activity
%and condition-dependent ROC curves.
%Set excludeSuppressed to 1 to exclude channels with stimulus-evoked
%suppression, i.e. blanco 13, 24, 42 and jack 49.
plotChFigs=1;
if useISI==1
    analysisType='ROC';
else
    analysisType='ROC_zero_one';
end
sglroc3IndividualChs=0;%set to 0 to read ROC values for individual channels and calculate mean ROC across channels; set to 1 to calculate ROCs based on pooled activity across channels
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
plotDiffC50_30=1;
calculateTangent=1;
if nargin<3||isempty(animals)
    animals=[{'blanco'} {'jack'}];
    % animals={'blanco'};
end
if nargin<4||isempty(areas)
    areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'}];
    areas=[{'v4_1'} {'v1_1'} {'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
end
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        channels=main_channels(animal,area);
        sessionNums=main_raw_sessions_final(animal,area,[],0);
        if exampleFig==1
            channels=51;
%             sessionNums=307;
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
            colmapText=colormap(jet(size(testContrast,2)));
            colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
            all_rocvals_sglroc3=[];
            all_rocvals=[];
            slopeNeuroNew=[];PNENew=[];diffPNENew=[];minRateNew=[];maxRateNew=[];chSSENew=[];
            slopeNeuroOld=[];PNEOld=[];diffPNEOld=[];minRateOld=[];maxRateOld=[];chSSEOld=[];
%             if useISI==1
                threshold82higher=[];
                threshold82lower=[];
%             end
            chPROBMAT={[]};
            chAUROC={[]};
            includeChCount=0;
            for chInd=1:length(channels)
                includeCh=1;
                if excludeSuppressed==1
                    if strcmp(animal,'blanco')&&channels(chInd)==13||channels(chInd)==24||channels(chInd)==42
                        includeCh=0;
                    elseif strcmp(animal,'jack')&&channels(chInd)==49
                        includeCh=0;
                    end
                end
                if includeCh==1
                    includeChCount=includeChCount+1;
                    if plotChFigs==1
                        figROC(includeChCount)=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
                        set(figROC(includeChCount), 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
                    end
                    sessionCounter=1;
                    for i=1:length(sessionNums)
                        rocvals=zeros(1,length(testContrast));
                        matFolder=['F:\PL\spikeData\',animal];
                        ROCmatChs=[];
                        chStr=[num2str(channels(chInd)),'_',num2str(sessionNums(i)),'_',num2str(sampleContrast),'.mat'];
                        matPath=fullfile(matFolder,chStr);
                        matExists=0;
                        if exist(matPath,'file')
                            matExists=1;
                        end
                        includeRows=includeSessionsAll(find(includeSessionsAll(:,1)==channels(chInd)),2);%include this session in analysis
                        includeRow=find(includeRows==sessionNums(i));
                        if matExists==1&&~isempty(includeRow)
                            valsText=['load ',matPath,' matarray'];
                            eval(valsText);
                            for condInd=1:length(testContrast)
                                if useISI==0
                                    if size(matarray{condInd,2},1)~=size(matarray{condInd,4},1)
                                        pause%if number of trials are not equal
                                    end
                                elseif useISI==1
                                    if size(matarray{condInd,3},1)~=size(matarray{condInd,4},1)
                                        pause%if number of trials are not equal
                                    end
                                end
                                if normalize==1
                                    maxAll=[];%find the highest firing rate across all conditions and trials, across both the sample and test presentation periods
                                    for n=1:size(matarray{condInd,4})
                                        if useISI==0
                                            maxAll=[maxAll length(matarray{condInd,2}{n})*1000/512 length(matarray{condInd,4}{n})*1000/512];
                                        elseif useISI==1
                                            temp3=matarray{condInd,3}{n}>512*2-256;%activity during second half of ISI
                                            spikes=matarray{condInd,3}{n}(temp3);
                                            temp3=spikes<512*2;
                                            spikes=spikes(temp3);
                                            maxAll=[maxAll length(spikes)*1000/256 length(matarray{condInd,4}{n})*1000/512];
                                        end
                                    end
                                    maxval=max(maxAll)/100;
                                else
                                    maxval=1;
                                end
                                chTestHigher=0;%tally number of trials with higher/lower test act for separate channels (not summing up activity across channels)
                                chTestLower=0;
                                actList2=[];
                                actList3=[];
                                actList4=[];
                                includeTrialCount=0;
                                for n=1:size(matarray{condInd,4})
                                    if useISI==0
                                        if (length(matarray{condInd,2}{n})*1000/512)/maxval~=(length(matarray{condInd,4}{n})*1000/512)/maxval%if test and sample act are not identical
                                            includeTrialCount=includeTrialCount+1;
                                            actList2(includeTrialCount)=(length(matarray{condInd,2}{n})*1000/512)/maxval;
                                            actList4(includeTrialCount)=(length(matarray{condInd,4}{n})*1000/512)/maxval;
                                            if actList2(includeTrialCount)<actList4(includeTrialCount)
                                                chTestHigher=chTestHigher+1;
                                            elseif actList2(includeTrialCount)>actList4(includeTrialCount)
                                                chTestLower=chTestLower+1;
                                            end
                                        end
                                    elseif useISI==1
                                        temp3=matarray{condInd,3}{n}>512*2-256;%activity during ISI
                                        spikes=matarray{condInd,3}{n}(temp3);
                                        temp3=spikes<512*2;
                                        spikes=spikes(temp3);
                                        if (length(spikes)/256*1000)/maxval~=(length(matarray{condInd,4}{n})*1000/512)/maxval
                                            includeTrialCount=includeTrialCount+1;
                                            actList3(includeTrialCount)=(length(spikes)/256*1000)/maxval;%find rate during second half of ISI
                                            actList4(includeTrialCount)=(length(matarray{condInd,4}{n})*1000/512)/maxval;
                                            if actList2(includeTrialCount)<actList4(includeTrialCount)
                                                chTestHigher=chTestHigher+1;
                                            elseif actList2(includeTrialCount)>actList4(includeTrialCount)
                                                chTestLower=chTestLower+1;
                                            end
                                        end
                                    end
                                end
                                %calculate individual channel PROBMAT:
                                chPROBMAT{includeChCount}(i,condInd)=chTestHigher/(chTestHigher+chTestLower);%sessions in rows, conditions in columns
                                %calculate individual channel ROC:
                                sessionsList=[];
                                chAUROC{chInd}(i,condInd)=sglroc3(actList4,actList2);%old method
                            end
                            if plotChFigs==1
                                %if ~isnan(PROBMATrocvals)
                                figure(figROC(includeChCount));
                                subplot(ceil(length(sessionNums)/5),5,sessionCounter);
                                xlim([0 max(testContrast)+10]);
                                if useISI==0
                                    plot(testContrast,chAUROC{includeChCount}(i,:),'bo');hold on%AUROC method
                                    [slopeNeuroOld,PNEOld,diffPNEOld,minRateOld,maxRateOld,chSSEOld]=weibull_fitting(chAUROC{includeChCount}(i,:),sampleContrast,testContrast,'old',sessionCounter,slopeNeuroOld,chSSEOld,PNEOld,minRateOld,maxRateOld,diffPNEOld,plotDiffC50_30,calculateTangent,useISI);
                                    plot(testContrast,chPROBMAT{includeChCount}(i,:),'ro');hold on%PROBMAT method
                                    if includeChCount==12&&i==16
                                        adjustCh=1;
                                    end
                                    [slopeNeuroNew,PNENew,diffPNENew,minRateNew,maxRateNew,chSSENew,threshold82higher,threshold82lower]=weibull_fitting(chPROBMAT{includeChCount}(i,:),sampleContrast,testContrast,'new',sessionCounter,slopeNeuroNew,chSSENew,PNENew,minRateNew,maxRateNew,diffPNENew,plotDiffC50_30,calculateTangent,useISI,threshold82higher,threshold82lower);
                                end
                                if sessionCounter==1
                                    xlabel('contrast (%)');
                                    ylabel('AUROC/PROBMAT');
                                end
                                sessionCounter=sessionCounter+1;
                                title(num2str(i),'FontSize',16);
                                %end
                            end
                        end
                    end
                    if useISI==0
                        subFolder='new_vs_old_individualchannels';
                    elseif useISI==1
                        subFolder='new_ROC_useISI_meanchannels';
                    end
                    if excludeSuppressed
                        subFolder=[subFolder,'_excludeSuppressed'];
                    end
                    if normalize
                        subFolder=[subFolder,'_normalised'];
                    end
                    if useISI==0
                        imagename=['ch_',num2str(channels(chInd)),'_ROCs_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
                    elseif useISI==1
                        imagename=['ch_',num2str(channels(chInd)),'_ROCs_old_new_useISI_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
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
                if useISI==0
                    [hS,pS,ciS,statsS]=ttest(slopeNeuroNew,slopeNeuroOld)
                    [hP,pP,ciP,statsP]=ttest(PNENew,PNEOld)
                    [hmin,pmin,cimin,statsmin]=ttest(minRateNew,minRateOld)
                    [hmax,pmax,cimax,statsmax]=ttest(maxRateNew,maxRateOld)
                end
                if useISI==0
                    matname=['ch_',num2str(channels(chInd)),'_ROCs_old_new_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
                elseif useISI==1
                    matname=['ch_',num2str(channels(chInd)),'_ROCs_old_new_useISI_',area,'_',num2str(sampleContrast),'_cutoff',num2str(cutoff*10)];
                end
                pathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder,matname);
                if useISI==0
                    saveText=['save ',pathname,'.mat all_rocvals_sglroc3 all_rocvals slopeNeuroNew PNENew diffPNENew minRateNew maxRateNew chSSENew slopeNeuroOld PNEOld diffPNEOld minRateOld maxRateOld chSSEOld hS pS ciS statsS hmin pmin cimin statsmin hP pP ciP statsP hmax pmax cimax statsmax threshold82higher threshold82lower'];
                elseif useISI==1
                    saveText=['save ',pathname,'.mat all_rocvals_sglroc3 all_rocvals slopeNeuroNew PNENew diffPNENew minRateNew maxRateNew chSSENew threshold82higher threshold82lower'];
                end
                eval(saveText);
            end
        end
    end
end
