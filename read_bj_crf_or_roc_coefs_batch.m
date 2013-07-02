function read_bj_crf_or_roc_coefs_batch(animal,area,analysisType,excludeSessHighSSE,excludeOutliers,comparisonType,plotLeastSquares,modifyMinMax,useISI,calcPartial)

%Written by Xing 13/06/13
%Just read correlation coefficients and calculate stats.
%excludeSessHighSSE: set to 1 to exclude sessions with poor Weibull fit
%excludeOutliers: set to 1 to exclude outlying data points
%modifyMinMax: set to 1 to read old minRate & maxRate values and calculate
%correct ones, set to 0 if values in file are already orrect.
%Set useISI=1 to read AUROC values that are calculated based on comparison of
%test to pre-test (rather than on test to sample) activity.
%Set calcPartial to 1 to run partial correlation instead of regular
%correlation analyses.
readData=0;
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
if nargin<7||isempty(plotLeastSquares)
    plotLeastSquares=[];
end
example_ch_54=0;
plotDiffC50_30=1;
writeCoefs=1;
slSigmaMultiple=[];
c50SigmaMultiple=[];
excludeSessions=[26 50 306 312 316 322:328 342 398 451];
test_epochs={0 512 512*2 512*3};durSpon=150;
channels=main_channels(animal,area);
[sampleContrasts testContrasts]=area_metadata(area);
psychoname=['psycho_constants_',area];
psychoPathname=fullfile('F:','PL','psycho_data',animal,psychoname);
if readData==1
    for i=1:length(channels)
        for sampleContrastsInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastsInd);
            testContrast=testContrasts(sampleContrastsInd,:);
            appendText=['_',num2str(sampleContrast)];
            for epoch=1:size(test_epochs,2)
                if strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC')&&epoch==4||strcmp(analysisType,'NVP')&&epoch==4||strcmp(analysisType,'ROC_diff')&&epoch==4||strcmp(analysisType,'ROC_zero_one')&&epoch==4||strcmp(analysisType,'NVP_zero_one')&&epoch==4
                    if epoch==1
                        periods=[-durSpon 0];
                    else
                        periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
                    end
                    for subPeriod=1:length(periods)-1
                        startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                        %read correlation coefficients, session numbers, and slope values
                        %from file:                        
                        chText=num2str(channels(i));
                        if excludeSessHighSSE==0
                            coefMatname=[chText,appendText,startEndTime,'_',analysisType,'_coefs_',area];
                        elseif excludeSessHighSSE==1
                            if excludeOutliers==0
                                coefMatname=[chText,appendText,startEndTime,'_',analysisType,'_coefs_',area,'_goodSSE'];
                            elseif excludeOutliers==1
                                coefMatname=[chText,appendText,startEndTime,'_',analysisType,'_coefs_',area,'_goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
                            end
                        end
                        if strcmp(analysisType,'CRF')
                            %                     if calcPartial==0
                            subFolderName=[analysisType,'_coef_mat'];
                            %                     elseif calcPartial==1
                            %                         subFolderName=[analysisType,'_coef_mat_partialcorr'];
                            %                     end
                            coefMatFolder=fullfile('F:','PL',analysisType,animal,subFolderName);
                            coefMatPathname=fullfile(coefMatFolder,coefMatname);
                        elseif strcmp(analysisType,'ROC_zero_one')
                            if calcPartial==1
                                coefMatPathname=['F:\PL\ROC_zero_one\partial_corr_4_factors\',animal,'\ROC_zero_one_coef_mat_partialcorr_4factors\',coefMatname];
                            else
                                coefMatPathname=['F:\PL\ROC_zero_one\',animal,'\ROC_zero_one_coef_mat\',coefMatname];
                            end
                        end
                        loadText=['load ',coefMatPathname,'.mat'];
                        eval(loadText)
                        plot_neurometric_coefs(animal,area,channels(i),appendText,startEndTime,slopeNeuro,c50,plotDiffC50_30,diffc50,minRate,maxRate,sessionSorted1,analysisType,example_ch_54,excludeSessHighSSE,excludeOutliers,writeCoefs,slSigmaMultiple,c50SigmaMultiple,calcPartial)
                    end
                end
            end
        end
    end
end
for sampleContrastsInd=1:length(sampleContrasts)
    sampleContrast=sampleContrasts(sampleContrastsInd);
    testContrast=testContrasts(sampleContrastsInd,:);
    appendText=['_',num2str(sampleContrast)];
    chCoefficients=[];
    chps=[];%p-vals
    includedChannels=[];
    c50s={[]};%the number of sessions incuded for each channel varies- as long as at least 80% of sessions have good slopes and PNE, channel is included
    slopes={[]};
    mins={[]};
    maxs={[]};
    chRowCouht=0;
    for i=1:length(channels)
        for epoch=4:4
            if strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC')&&epoch==4||strcmp(analysisType,'NVP')&&epoch==4||strcmp(analysisType,'ROC_diff')&&epoch==4||strcmp(analysisType,'ROC_zero_one')&&epoch==4||strcmp(analysisType,'NVP_zero_one')&&epoch==4
                if epoch==1
                    periods=[-durSpon 0];
                else
                    periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
                end
                for subPeriod=1:length(periods)-1
                    startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                    %read correlation coefficients, session numbers, and slope values
                    %from file:
                    chText=num2str(channels(i));
                    if excludeSessHighSSE==0
                        coefMatname=[chText,appendText,startEndTime,'_',analysisType,'_coefs_',area,'.mat'];
                    elseif excludeSessHighSSE==1
                        if excludeOutliers==0
                            coefMatname=[chText,appendText,startEndTime,'_',analysisType,'_coefs_',area,'_goodSSE','.mat'];
                        elseif excludeOutliers==1
                            coefMatname=[chText,appendText,startEndTime,'_',analysisType,'_coefs_',area,'_goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple),'.mat'];
                        end
                    end
                    if calcPartial==0
                        subFolderName=[analysisType,'_coef_mat'];
                    elseif calcPartial==1
                        subFolderName=[analysisType,'_coef_mat_partialcorr'];
                    end
                    coefMatFolder=fullfile('F:','PL',analysisType,animal,subFolderName);
                    coefMatPathname=fullfile(coefMatFolder,coefMatname);
                    if exist(coefMatPathname,'file')
                        loadText=['load ',coefMatPathname,' coefficients c50 slopeNeuro minRate maxRate sessionSorted1'];
                        eval(loadText)
                        chRowCouht=chRowCouht+1;
                        chCoefficients=[chCoefficients;coefficients(1,:)];
                        chps=[chps;coefficients(2,:)];
                        includedChannels=[includedChannels channels(i)];
                        c50s{chRowCouht}=c50;%record C50 values, will need to tally number of channels where C50 shifts towards or away from 30%
                        slopes{chRowCouht}=slopeNeuro;%record slope values, will need to tally number of channels where slope becomes steeper or shallower
                        mins{chRowCouht}=minRate;%record slope values, will need to tally number of channels where slope becomes steeper or shallower
                        maxs{chRowCouht}=maxRate;%record slope values, will need to tally number of channels where slope becomes steeper or shallower
                    end
                end
            end
        end
    end
    chCoefficients;
    sigDirection=[];
    sigChDirSlope=cell(3,1);
    dir30=[0;0];%towards or away from 30
    dirSlope=[0;0;0];%steeper or shallower
    dirMin=[0;0];%increase or decrease in min
    dirMax=[0;0];%increase or decrease in max
    sigChNames=cell(5,1);
    sigChNamesDir=cell(10,1);
    sigChDir30=cell(2,1);
    sigChDirMin=cell(2,1);
    sigChDirMax=cell(2,1);
    for paramInd=1:size(chCoefficients,2)
%         if calcPartial==0
%             sigInd=find(chps(:,paramInd)<0.05/size(chps,1));
%         elseif calcPartial==1
            sigInd=find(chps(:,paramInd)<0.05);
%         end
        sigChNames(paramInd,1)={includedChannels(sigInd)};%compile list of channel numbers
        sigChs(paramInd,1)=length(sigInd);
        sigCoefs=chCoefficients(sigInd,paramInd);
        sigDirection=[sigDirection;sum(sigCoefs>0);sum(sigCoefs<0)];%number of positive & negative changes
        sigChNamesDir(paramInd*2-1)={sigChNames{paramInd,1}(sigCoefs>0)};%list of channel numbers w positive 
        sigChNamesDir(paramInd*2)={sigChNames{paramInd,1}(sigCoefs<0)};%& negative changes
        if ~isempty(sigInd)
            if paramInd==1%sig change in slope
                for chRowInd=1:length(sigInd)
                    sigSlopes=slopes{sigInd(chRowInd)};
                    lengthThird=floor(size(sigSlopes,2)/3);
                    if abs(mean(sigSlopes))>0&&mean(sigSlopes(1:lengthThird))<mean(sigSlopes(end-lengthThird+1:end))||abs(mean(sigSlopes))<0&&mean(sigSlopes(1:lengthThird))>mean(sigSlopes(end-lengthThird+1:end))
                        dirSlope(1,1)=dirSlope(1)+1;%steeper
                        sigChDirSlope{1,1}=[sigChDirSlope{1,1} sigChNames{paramInd,1}(chRowInd)];
                    elseif abs(mean(sigSlopes))>0&&mean(sigSlopes(1:lengthThird))>mean(sigSlopes(end-lengthThird+1:end))||abs(mean(sigSlopes))<0&&mean(sigSlopes(1:lengthThird))<mean(sigSlopes(end-lengthThird+1:end))
                        dirSlope(2,1)=dirSlope(2)+1;%shallower
                        sigChDirSlope{2,1}=[sigChDirSlope{2,1} sigChNames{paramInd,1}(chRowInd)];
                    elseif mean(sigSlopes(1:lengthThird))<0&&mean(sigSlopes(end-lengthThird+1:end))>0||mean(sigSlopes(1:lengthThird))>0&&mean(sigSlopes(end-lengthThird+1:end))<0
                        dirSlope(3,1)=dirSlope(3)+1;
                        sigChDirSlope{3,1}=[sigChDirSlope{3,1} sigChNames{paramInd,1}(chRowInd)];%switch from one polarity to another (i.e. +ve to -ve or -ve to +ve)
                    end
                end
            elseif paramInd==2%sig change in C50
                for chRowInd=1:length(sigInd)
                    sigchC50s=c50s{sigInd};
                    lengthThird=floor(size(sigchC50s,2)/3);
                    if mean(sigchC50s(1:lengthThird))<30&&mean(sigchC50s(1:lengthThird))<mean(sigchC50s(end-lengthThird+1:end))||mean(sigchC50s(1:lengthThird))>30&&mean(sigchC50s(1:lengthThird))>mean(sigchC50s(end-lengthThird+1:end))
                        dir30(1,1)=dir30(1)+1;%towards 30%
                        sigChDir30{1,1}=[sigChDir30{1,1} sigChNames{paramInd,1}(chRowInd)];
                    end
                    if mean(sigchC50s(1:lengthThird))<30&&mean(sigchC50s(1:lengthThird))>mean(sigchC50s(end-lengthThird+1:end))||mean(sigchC50s(1:lengthThird))>30&&mean(sigchC50s(1:lengthThird))<mean(sigchC50s(end-lengthThird+1:end))
                        dir30(2,1)=dir30(2)+1;%away from 30%
                        sigChDir30{2,1}=[sigChDir30{2,1} sigChNames{paramInd,1}(chRowInd)];
                    end
                end
            elseif paramInd==3%sig change in min
                for chRowInd=1:length(sigInd)
                    sigchMins=mins{sigInd};
                    lengthThird=floor(size(sigchMins,2)/3);
                    if mean(sigchMins(1:lengthThird))<mean(sigchMins(end-lengthThird+1:end))
                        dirMin(1,1)=dirMin(1)+1;
                        sigChDirMin{1,1}=[sigChDirMin{1,1} sigChNames{paramInd,1}(chRowInd)];%increase
                    end
                    if mean(sigchMins(1:lengthThird))>mean(sigchMins(end-lengthThird+1:end))
                        dirMin(2,1)=dirMin(2)+1;
                        sigChDirMin{2,1}=[sigChDirMin{2,1} sigChNames{paramInd,1}(chRowInd)];%decrease
                    end
                end
            elseif paramInd==4%sig change in max
                for chRowInd=1:length(sigInd)
                    sigchMaxs=maxs{sigInd};
                    lengthThird=floor(size(sigchMaxs,2)/3);
                    if mean(sigchMaxs(1:lengthThird))<mean(sigchMaxs(end-lengthThird+1:end))
                        dirMax(1,1)=dirMax(1)+1;
                        sigChDirMax{1,1}=[sigChDirMax{1,1} sigChNames{paramInd,1}(chRowInd)];%increase
                    end
                    if mean(sigchMaxs(1:lengthThird))>mean(sigchMaxs(end-lengthThird+1:end))
                        dirMax(2,1)=dirMax(2)+1;
                        sigChDirMax{2,1}=[sigChDirMax{2,1} sigChNames{paramInd,1}(chRowInd)];%decrease
                    end
                end
            end
        end
    end
    sigChs%1 slope 2 C50 3 min 4 max 5 neuro vs psycho slopes
    sigDirection
    dir30
    formattedTable=[sigDirection(1:2);dirSlope;sigDirection(3:4);dir30;dirMin;dirMax;sigDirection(end-1:end)];
    formattedChsTable=[sigChNamesDir(1:2);sigChDirSlope;sigChNamesDir(3:4);sigChDir30;sigChDirMin;sigChDirMax;sigChNamesDir(end-1:end)];
    tallyMatname=['tallySigCoefs_',area,appendText];
    if strcmp(analysisType,'CRF')
        %                     if calcPartial==0
        subFolderName=[analysisType,'_coef_mat'];
        %                     elseif calcPartial==1
        %                         subFolderName=[analysisType,'_coef_mat_partialcorr'];
        %                     end
        coefMatFolder=fullfile('F:','PL',analysisType,animal,subFolderName);
        coefMatPathname=fullfile(coefMatFolder,tallyMatname);
    elseif strcmp(analysisType,'ROC_zero_one')
        if calcPartial==0
            coefMatPathname=['F:\PL\ROC_zero_one\',animal,'\ROC_zero_one_coef_mat\',tallyMatname];
        elseif calcPartial==1
            coefMatPathname=['F:\PL\ROC_zero_one\partial_corr_4_factors\',animal,'\ROC_zero_one_coef_mat_partialcorr_4factors\',tallyMatname];
        end
    end
    saveText=['save ',coefMatPathname,'.mat sigChs sigDirection dir30 sigCoefs sigChNamesDir sigChDir30 formattedTable formattedChsTable'];
    eval(saveText)
end

