function read_bj_crf_or_roc_batch(animal,area,analysisType,excludeSessHighSSE,excludeOutliers)

%Written by Xing 06/03/13
%to calculate and write correlation coefficients
%excludeSessHighSSE: set to 1 to exclude sessions with poor Weibull fit
%excludeOutliers: set to 1 to exclude outlying data points 
if strcmp(analysisType,'ROC_diff')
    comparisonType=1;
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
excludeSessions=[26 50 306 312 316 322:328 342];
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
            if strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC')&&epoch==4||strcmp(analysisType,'NVP')&&epoch==4||strcmp(analysisType,'ROC_diff')&&epoch==4
                if epoch==1
                    periods=[-durSpon 0];
                else
                    periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
                end
                for subPeriod=1:length(periods)-1
                    startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                    if strcmp(analysisType,'ROC_diff')%load ROC values from previous sglroc3 method
                        matName=['ROC','_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        if strcmp(area,'v4_1')||strcmp(area,'v1_1')
                            if strcmp(ROC1,'ROC_sglroc3_correct_only')
                                matName=['ROC','_Ch',num2str(channels(i)),'_',num2str(sampleContrast),'_1058_to_1587.mat'];
                            end
                        end
                        matPath=fullfile('F:','PL',ROC1,animal,area,matName);
                        loadText=['load ',matPath,' ','ROCmat'];
                        eval(loadText);
                        dataArray=ROCmat;
                        includeMatch=[];
                        for includeInd=1:size(dataArray,1)
                            if sum(excludeSessions==cell2mat(dataArray(includeInd,1)))==0
                                includeMatch=[includeMatch includeInd];
                            end
                        end
                        dataArray=dataArray(includeMatch,:);
                        read_bj_crf_or_roc(dataArray,channels(i),psychoPathname,testContrast,sampleContrast,animal,area,startEndTime,'ROC_diff2',excludeSessHighSSE,excludeOutliers,rootFolder)
                    end
                    if strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC')
                        matName=[analysisType,'_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        matPath=fullfile('F:','PL',analysisType,animal,area,matName);
                        loadText=['load ',matPath,' ',analysisType,'mat'];
                    elseif strcmp(analysisType,'ROC_diff')
                        matName=['ROC','_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        matPath=fullfile('F:','PL',ROC2,animal,area,matName);
                        loadText=['load ',matPath,' ','ROCmat'];
                    elseif strcmp(analysisType,'NVP')
                        matName=['ROC_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        matPath=fullfile('F:','PL','ROC',animal,area,matName);
                        loadText=['load ',matPath,' ROCmat'];
                    end
                    eval(loadText);
                    if strcmp(analysisType,'CRF')
                        dataArray=CRFmat;
                    elseif strcmp(analysisType,'ROC')
                        dataArray=ROCmat;
                    elseif strcmp(analysisType,'ROC_diff')
                        dataArray=ROCmat;
                    elseif strcmp(analysisType,'NVP')
                        dataArray=ROCmat;
                    end
                    includeMatch=[];
                    for includeInd=1:size(dataArray,1)
                        if sum(excludeSessions==cell2mat(dataArray(includeInd,1)))==0
                            includeMatch=[includeMatch includeInd];
                        end
                    end
                    dataArray=dataArray(includeMatch,:);
                    read_bj_crf_or_roc(dataArray,channels(i),psychoPathname,testContrast,sampleContrast,animal,area,startEndTime,analysisType,excludeSessHighSSE,excludeOutliers,rootFolder)
                end
            end
        end
    end
end
