function add_bj_crf_or_roc_batch_v4_0(animal,area,analysisType)

%Written by Xing 13/05/14
%Sort sessions and prepend data from V4_0_1 and V4_0_2 (with 8 and 12
%conditions, respectively) to ROCmat.
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
if nargin<7||isempty(plotLeastSquares)
    plotLeastSquares=[];
end
test_epochs={0 512 512*2 512*3};durSpon=150;
channels=main_channels(animal,area);
[sampleContrasts testContrasts]=area_metadata(area);
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
                    if strcmp(analysisType,'ROC')
                        matName=[analysisType,'_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        matPath=fullfile('F:','PL',analysisType,animal,area,matName);
                        loadText=['load ',matPath,' ',analysisType,'mat'];
                        saveText=['save ',matPath,' ',analysisType,'mat'];
                        eval(loadText);
                        [sessionsSorted ind]=sort(cell2mat(ROCmat(:,1)));
                        ROCmat=ROCmat(ind,:);
                        ROCmatFull=ROCmat;
                        matPath=fullfile('F:','PL',analysisType,animal,'v4_0_2',matName);
                        loadText=['load ',matPath,' ',analysisType,'mat'];
                        eval(loadText);
                        ROCmatFull=[ROCmat;ROCmatFull];
                        matPath=fullfile('F:','PL',analysisType,animal,'v4_0_1',matName);
                        loadText=['load ',matPath,' ',analysisType,'mat'];
                        eval(loadText);
                        ROCmatFull=[ROCmat;ROCmatFull];
                        ROCmat=ROCmatFull;
                        [sessionsSorted ind]=sort(cell2mat(ROCmat(:,1)));
                        ROCmat=ROCmat(ind,:);
                        %check against Mehdi's list of included sessions
                        loadMehdiSessions=['load F:\PL\mehdi_snr\',animal,'_SNR_V4_N.mat'];
                        eval(loadMehdiSessions)
                        if length(cell2mat(ROCmat(:,1)))==length(REC_PEN)-1
                            if sum(cell2mat(ROCmat(:,1))'==REC_PEN(1:end-1))==size(ROCmat,1)
                                eval(saveText)
                            end
                        else
                            pause
                        end
                    elseif strcmp(analysisType,'CRF')
                        matName=[analysisType,'_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                        matPath=fullfile('F:','PL',analysisType,animal,area,matName);
                        loadText=['load ',matPath,' ',analysisType,'mat'];
                        saveText=['save ',matPath,' ',analysisType,'mat CRFmat'];
                        eval(loadText);
                        [sessionsSorted ind]=sort(cell2mat(CRFmat(:,1)));
                        CRFmat=CRFmat(ind,:);
                        CRFmatFull=CRFmat;
                        matPath=fullfile('F:','PL',analysisType,animal,'v4_0_2',matName);
                        loadText=['load ',matPath,' ',analysisType,'mat'];
                        eval(loadText);
                        CRFmatFull=[CRFmat;CRFmatFull];
                        matPath=fullfile('F:','PL',analysisType,animal,'v4_0_1',matName);
                        loadText=['load ',matPath,' ',analysisType,'mat'];
                        eval(loadText);
                        CRFmatFull=[CRFmat;CRFmatFull];
                        CRFmat=CRFmatFull;
                        [sessionsSorted ind]=sort(cell2mat(CRFmat(:,1)));
                        CRFmat=CRFmat(ind,:);
                        %check against Mehdi's list of included sessions
                        loadMehdiSessions=['load F:\PL\mehdi_snr\',animal,'_SNR_V4_N.mat'];
                        eval(loadMehdiSessions)
                        if length(cell2mat(CRFmat(:,1)))==length(REC_PEN)-1
                            if sum(cell2mat(CRFmat(:,1))'==REC_PEN(1:end-1))==size(CRFmat,1)
                                eval(saveText)
                            end
                        else
                            pause
                        end               
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
                end
            end
        end
    end
end
