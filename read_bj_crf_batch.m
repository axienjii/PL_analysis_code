function read_bj_crf_batch(animal,area)

%Written by Xing 06/03/13
%to calculate and write correlation coefficients
test_epochs={0 529 529*2 529*3};durSpon=150;
channels=main_channels(animal,area);
[sampleContrasts testContrasts]=area_metadata(area);
psychoname=['psycho_constants_',area];
psychoPathname=fullfile('F:','PL','psycho_data',animal,psychoname);
for i=1:length(channels)
    for sampleContrastsInd=1:length(sampleContrasts)
        sampleContrast=sampleContrasts(sampleContrastsInd);
        testContrast=testContrasts(sampleContrastsInd,:);    
        for epoch=1:size(test_epochs,2)
            if epoch==1
                periods=[-durSpon 0];
            else
                periods=[test_epochs{epoch-1} test_epochs{epoch}(1)];
            end
            for subPeriod=1:length(periods)-1
                startEndTime=['_',num2str(periods(subPeriod)),'_to_',num2str(periods(subPeriod+1))];
                CRFmatName=['CRF_Ch',num2str(channels(i)),'_',num2str(sampleContrast),startEndTime,'.mat'];
                CRFmatPath=fullfile('F:','PL','CRF',animal,area,CRFmatName);
                loadText=['load ',CRFmatPath,' CRFmat'];
                eval(loadText);
                read_bj_crf(CRFmat,channels(i),psychoPathname,testContrast,sampleContrast,animal,area,startEndTime)
            end
        end
    end
end
