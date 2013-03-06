function read_bj_crf_batch(animal,area)

%Written by Xing 06/03/13
%to calculate and write correlation coefficients:
channels=main_channels(animal,area);
[sampleContrasts testContrasts]=area_metadata(area);
psychoname=['psycho_constants_',area];
psychoPathname=fullfile('F:','PL','psycho_data',animal,psychoname);
for i=1:length(channels)
    for sampleContrastsInd=1:length(sampleContrasts)
        sampleContrast=sampleContrasts(sampleContrastsInd);
        testContrast=testContrasts(sampleContrastsInd,:);
            CRFmatName=['CRF_Ch',num2str(channels(i)),'_',num2str(sampleContrast),'.mat'];
        CRFmatPath=fullfile('F:','PL','CRF',animal,area,CRFmatName);
        loadText=['load ',CRFmatPath,' CRFmat'];
        eval(loadText);
        read_bj_crf(CRFmat,channels(i),psychoPathname,testContrast,sampleContrast,animal,area)
    end
end
