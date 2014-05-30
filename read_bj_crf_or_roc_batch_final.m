function read_bj_crf_or_roc_batch_final(animal,area,analysisType,excludeSessHighSSE,excludeOutliers,comparisonType,plotLeastSquares,modifyMinMax,useISI,calcPartial)

%Written by Xing 06/03/13
%to calculate and write correlation coefficients
%excludeSessHighSSE: set to 1 to exclude sessions with poor Weibull fit
%excludeOutliers: set to 1 to exclude outlying data points 
%modifyMinMax: set to 1 to read old minRate & maxRate values and calculate
%correct ones, set to 0 if values in file are already orrect.
%Set useISI=1 to read AUROC values that are calculated based on comparison of
%test to pre-test (rather than on test to sample) activity.
%Set calcPartial to 1 to run partial correlation instead of regular
%correlation analyses.
if strcmp(analysisType,'ROC_diff')
    switch(comparisonType)
        case(1)%compare ROC values between two methods, using all trials
            ROC1='ROC_sglroc3';
            ROC2='ROC';
        case(2)%compare ROC values between two methods, using only correct trials
            ROC1='ROC_sglroc3_correct_only';
            ROC2='ROC_correct_only';
        case(3)%compare ROC values between correct versus all trials, using new method
            ROC1='ROC_correct_only';
            ROC2='ROC';
    end
end
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
if nargin<7||isempty(plotLeastSquares)
    plotLeastSquares=[];
end
excludeSessions=[26 50 306 312 316 322:328 342 398 451];%before blu-ray data added, from sessions for which raw data is missing
excludeSessions=[26 50 312 342 398 451];
test_epochs={0 512 512*2 512*3};durSpon=150;
channels=main_channels(animal,area);
[sampleContrasts testContrasts]=area_metadata(area);
psychoname=['psycho_constants_',area];
psychoPathname=fullfile('F:','PL','psycho_data',animal,psychoname);
for i=1:length(channels)
    for sampleContrastsInd=1:length(sampleContrasts)
        sampleContrast=sampleContrasts(sampleContrastsInd);
        testContrast=testContrasts(sampleContrastsInd,:);    
        for epoch=1:size(test_epochs,2)
            if strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC')&&epoch==4||strcmp(analysisType,'NVP')&&epoch==4||strcmp(analysisType,'ROC_diff')&&epoch==4||strcmp(analysisType,'ROC_zero_one')&&epoch==4||strcmp(analysisType,'NVP_zero_one')&&epoch==4
                if epoch==1
                    periods=[-durSpon 0];
                else
                    periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
                end
                for subPeriod=1:length(periods)-1
                    startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                    if strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC')
                        matName=[analysisType,'_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        matPath=fullfile('F:','PL',analysisType,animal,area,matName);
                        loadText=['load ',matPath,' ',analysisType,'mat'];
                        saveText=['save ',matPath,' ',analysisType,'mat ROCmat'];
                    elseif strcmp(analysisType,'ROC_zero_one')
                        matName=['ROC_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        matPath=fullfile('F:','PL',analysisType,'ROC',animal,area,matName);
                        loadText=['load ',matPath,' ','ROCmat'];
                    elseif strcmp(analysisType,'ROC_diff')
                        matName=['ROC','_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        matPath=fullfile('F:','PL',ROC2,animal,area,matName);
                        loadText=['load ',matPath,' ','ROCmat'];
                    elseif strcmp(analysisType,'NVP')
                        matName=['ROC_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        matPath=fullfile('F:','PL','ROC',animal,area,matName);
                        loadText=['load ',matPath,' ROCmat'];
                    elseif strcmp(analysisType,'NVP')||strcmp(analysisType,'NVP_zero_one')
                        matName=['ROC_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        matPath=fullfile('F:','PL','ROC_zero_one','ROC',animal,area,matName);
                        loadText=['load ',matPath,' ROCmat'];
                    end
                    eval(loadText);
                    if strcmp(analysisType,'CRF')
                        [sessionsSorted ind]=sort(cell2mat(CRFmat(:,1)));
                        CRFmat=CRFmat(ind,:);
                        eval(saveText)
                        dataArray=CRFmat;
                    elseif strcmp(analysisType,'ROC')||strcmp(analysisType,'ROC_zero_one')
                        [sessionsSorted ind]=sort(cell2mat(ROCmat(:,1)));
                        ROCmat=ROCmat(ind,:);
                        eval(saveText)
                        dataArray=ROCmat;
                    elseif strcmp(analysisType,'ROC_diff')
                        dataArray=ROCmat;
                    elseif strcmp(analysisType,'NVP')||strcmp(analysisType,'NVP_zero_one')
                        dataArray=ROCmat;
                    end
                    %basic session inclusion:
                    includeMatch=[];
                    for includeInd=1:size(dataArray,1)
                        if sum(excludeSessions==cell2mat(dataArray(includeInd,1)))==0
                            includeMatch=[includeMatch includeInd];
                        end
                    end
                    dataArray=dataArray(includeMatch,:);
                    %more refined session inclusion, based on Mehdi's SNR
                    %criteria:
                    loadText=['load F:\PL\SNR\',animal,'_',area,'_includeSessions.mat goodChannels goodSessions'];
                    eval(loadText)
                    chInd=find(goodChannels==channels(i));
                    if ~isempty(chInd)
                        includeSessions=goodSessions{chInd};
                        includeMatch=[];
                        for includeInd=1:size(dataArray,1)
                            if sum(includeSessions==cell2mat(dataArray(includeInd,1)))==1
                                includeMatch=[includeMatch includeInd];
                            end
                        end
                        dataArray=dataArray(includeMatch,:);
                        read_bj_crf_or_roc_final(dataArray,channels(i),psychoPathname,testContrast,sampleContrast,animal,area,startEndTime,analysisType,excludeSessHighSSE,excludeOutliers,rootFolder,plotLeastSquares,modifyMinMax,useISI,calcPartial)
                    end
                end
            end
        end
    end
end
