function read_bj_crf_or_roc_batch(animal,area,analysisType,excludeSessHighSSE,excludeOutliers)

%Written by Xing 06/03/13
%to calculate and write correlation coefficients
%excludeSessHighSSE: set to 1 to exclude sessions with poor Weibull fit
%excludeOutliers: set to 1 to exclude outlying data points 
excludeSessions=[26 50 306 312 316 322:328 342];
if strcmp(area,'v4_1')||strcmp(area,'v1_1')
    test_epochs={0 529 529*2 529*3};durSpon=150;
elseif strncmp(area,'v1_2',4)
    test_epochs={0 512 512*2 512*3};durSpon=150;
end
channels=main_channels(animal,area);
[sampleContrasts testContrasts]=area_metadata(area);
psychoname=['psycho_constants_',area];
psychoPathname=fullfile('F:','PL','psycho_data',animal,psychoname);
for i=1:length(channels)
    for sampleContrastsInd=1:length(sampleContrasts)
        sampleContrast=sampleContrasts(sampleContrastsInd);
        testContrast=testContrasts(sampleContrastsInd,:);    
        for epoch=1:size(test_epochs,2)
            if strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC')&&epoch==4||strcmp(analysisType,'NVP')&&epoch==4
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
                    read_bj_crf_or_roc(dataArray,channels(i),psychoPathname,testContrast,sampleContrast,animal,area,startEndTime,analysisType,excludeSessHighSSE,excludeOutliers)
                end
            end
        end
    end
end
