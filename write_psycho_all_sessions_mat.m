function write_psycho_all_sessions_mat
animals=[{'blanco'} {'jack'}];
areas=[{'v4_1'} {'v1_1'}];
areas=[{'v4_1'} {'v4_2'} {'v1_1'} {'v1_2'} {'v1_2_1'} {'v1_2_2'} {'v1_2_3'}];
test_epochs={0 512 512*2 512*3};durSpon=150;
durSpon=150;%length of period prior to sample onset from which spontaneous rates are calculated. Can take on a value of up to 512 ms.
minTrials=10;%set value of minumum number of trials for inclusion of session
subPeriod=1;
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        [sampleContrasts testContrasts]=area_metadata(area);
        for sampleContrastsInd=1:length(sampleContrasts)
            sampleContrast=sampleContrasts(sampleContrastsInd);
            testContrast=testContrasts(sampleContrastsInd,:);
            loadText=['load F:\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_',area,'_',num2str(sampleContrast),'.mat'];
            eval(loadText);
            psychoAll=[allMeanPerf(:,1) allMeanPerf(:,3:length(testContrast)+6)];
            saveText=['save F:\PL\psycho_data\',animal,'\',area,'_psycho_all_sessions.mat psychoAll'];
            eval(saveText)
        end
    end
end