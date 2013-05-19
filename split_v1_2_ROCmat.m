function split_v1_2_ROCmat(ROCtype)
animals=[{'blanco'} {'jack'}];
areas={'v1_2_1' 'v1_2_2' 'v1_2_3'};
if strcmp(ROCtype,'sglroc3')
    ROCtext='_sglroc3';
else
    ROCtext=[];
end
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        sessions=main_raw_sessions_final(animal,area,[],0);
        for sampleInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleInd);
            channels = main_channels(animal,area);
            for i=1:length(channels)
                loadText=['load F:\PL\ROC',ROCtext,'\',animal,'\v1_2\ROC_Ch',num2str(channels(i)),'_',num2str(sampleContrast),'_1024_to_1536'];
                eval(loadText);
                sessionsList=[];
                for j=1:size(ROCmat,1)
                    sessionsList=[sessionsList ROCmat{j,1}];
                end
                rowsInd=[];
                for j=1:length(sessions)
                    rowInd=find(sessionsList==sessions(j));
                    rowsInd=[rowsInd rowInd];
                end
                ROCmat=ROCmat(rowsInd,:);
                folder=['F:\PL\ROC',ROCtext,'\',animal,'\',area];
                if ~exist(folder,'dir')
                    mkdir(folder)
                end
                saveText=['save F:\PL\ROC',ROCtext,'\',animal,'\',area,'\ROC_Ch',num2str(channels(i)),'_',num2str(sampleContrast),'_1024_to_1536 ROCmat'];
                eval(saveText);
            end
        end
    end
end