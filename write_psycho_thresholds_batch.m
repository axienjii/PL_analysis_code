function write_psycho_thresholds_batch(roving,excludeSessHighSSE,excludeOutliers,analysisType,useISI)
plotLeastSquares=0;
onExternalHD=1;
if onExternalHD==0
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
% analysisType='psycho';
% analysisType='psycho_param';
startEndTime='wholetrial';
animalTexts=[{'subject B'} {'subject J'}];
animals=[{'blanco'} {'jack'}];
SSEcutoff=0.09;
manual_cutoff=100;
threshSigmaMultiple=3;
excludeSessions=[26 50 306 312 316 322:328 342];
test_epochs={0 512 512*2 512*3};
if roving==0
    areaTexts=[{'V4'} {'V1'}];
    areas=[{'v4_1'} {'v1_1'} {'v4_2'}];
elseif roving==1
    areaTexts={'V1 roving data'};
    areas={'v1_2_1' 'v1_2_2' 'v1_2_3'};
end
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        sessions=main_raw_sessions_final_psycho(animal,area,[],0);
        sessionSorted1=[];
        for sampleInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleInd);
            testContrast=testContrasts(sampleInd,:);
            appendText=['_',num2str(sampleContrast)];
            if roving==0
                loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\',area,'_psycho_all_sessions.mat psychoAll'];
            elseif roving==1
                loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\v1_2_psycho_all_sessions',appendText,'.mat psychoAll'];
            end
            eval(loadText)
            psychomat=[];
            sessionSorted1=[];
            for sessionInd=1:length(sessions)
                for rowInd=1:size(psychoAll,1)
                    if psychoAll(rowInd,1)==sessions(sessionInd)
                        psychomat=[psychomat;{psychoAll(rowInd,1)} {test_epochs} {psychoAll(rowInd,2:length(testContrast)+1)}];
                        sessionSorted1=[sessionSorted1 psychoAll(rowInd,1)];
                    end
                end
            end
            numsessions=length(sessionSorted1);
            datamat=psychomat;
            saveText=['save ',rootFolder,'\PL\psycho_data\',animal,'\',area,'_psycho_array.mat psychomat'];
            eval(saveText)
            chSSE=zeros(length(sessionSorted1),2);
            SSEMatFileName=[area,appendText,startEndTime,'_SSE'];
            SSEMatFolder=fullfile(rootFolder,'PL',analysisType,animal,'SSE_mat_files');
            if ~exist(SSEMatFolder,'dir')
                mkdir(SSEMatFolder);
            end
            SSEMatPath=fullfile(SSEMatFolder,SSEMatFileName);
            psychoThresholdMatName=[area,appendText,'wholetrial_psyThreshold'];
            psychoThresholdMatFolder=fullfile(rootFolder,'PL',analysisType,animal,'psyThreshold_mat');
            if ~exist(psychoThresholdMatFolder,'dir')
                mkdir(psychoThresholdMatFolder);
            end
            psychoThresholdMatPathname=fullfile(psychoThresholdMatFolder,psychoThresholdMatName);
            if excludeSessHighSSE==1
                loadText=['load ',SSEMatPath,' chSSE'];
                eval(loadText);
                datamat=psychomat;
                SSEcutoff=mean(chSSE(:,2))+std(chSSE(:,2));
                ind=(chSSE(:,2)<SSEcutoff);
                sessionSorted1=sessionSorted1(ind);
                datamat=datamat(ind,:);
                numsessions=length(sessionSorted1);
                if excludeOutliers==1
                    loadText=['load ',psychoThresholdMatPathname];
                    eval(loadText)
                    tHsigma=std(threshold82higher);
                    tHoutliers=abs((threshold82higher-mean(threshold82higher)))>threshSigmaMultiple*tHsigma;
                    thHOutliersHighcut=threshold82higher>manual_cutoff;%manual exclusion for obvious outliers
                    if strcmp(analysisType,'psycho_zero_one')||strcmp(analysisType,'psycho_param_zero_one')
                        tLsigma=std(threshold82lower);
                        tLoutliers=abs((threshold82lower-mean(threshold82lower)))>threshSigmaMultiple*tLsigma;
                        thLOutliersHighcut=threshold82lower>manual_cutoff;%manual exclusion for obvious outliers (above 100%)
                        ind=tLoutliers+tHoutliers+thLOutliersHighcut+thHOutliersHighcut;%find sessions where lower and/or higher contrast threshold values are outliers (union)
                    else
                        ind=tHoutliers+thHOutliersHighcut;%find sessions where lower and/or higher contrast threshold values are outliers (union)
                    end
                    if sum(ind)>0
                        sessionSorted1=sessionSorted2(~ind);%keep sessions that do not have outliers
                        datamat=datamat(~ind,:);
                        numsessions=length(sessionSorted1);
                        if strcmp(analysisType,'psycho_zero_one')||strcmp(analysisType,'psycho_param_zero_one')
                            threshold82lower=threshold82lower(~ind);
                        else
                            threshold82lower=[];
                        end
                        threshold82higher=threshold82higher(~ind);
                    end
                end
            end
            [slopeNeuro,c50,diffc50,minRate,maxRate,chSSE,yLimData,threshold82lower,threshold82higher]=plot_CRF_or_ROC_across_sessions(animal,area,analysisType,datamat,analysisType,numsessions,sessionSorted1,sampleContrast,testContrast,1,1,excludeSessHighSSE,excludeOutliers,SSEMatPath,startEndTime,psychoThresholdMatPathname,[],[],threshSigmaMultiple,rootFolder,plotLeastSquares,useISI);
        end
    end
end
close all