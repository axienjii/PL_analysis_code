function PL_LFP_copy_files_daniel
%Written by Xing 10/04/14. Copies LFP .mat files to a secondary location.
newLocation='J:\PL';
animals=[{'blanco'} {'jack'}];
areas=[{'v4_1'} {'v1_1'}];
for animalInd=1:length(animals)
    animal=animals{animalInd};
    for areaInd=1:length(areas)
        area=areas{areaInd};
        channels=main_channels(animal,area);
        sessions=main_raw_sessions_final(animal,area,[],0);
        for j=1:length(sessions)
            for i=1:length(channels)
                originalPath=['M:\Xing\pl_LFP\',animal,'\',num2str(sessions(j)),'\SGL_trial_LFP_-512_1536_ch',num2str(channels(i)),'.mat'];
                newFolder=[newLocation,'\',animal,'\',num2str(sessions(j))];
                if ~exist(newFolder,'dir')
                    mkdir(newFolder)
                end
                newPath=[newFolder,'\SGL_trial_LFP_-512_1536_ch',num2str(channels(i)),'.mat'];
                copyfile(originalPath,newPath);
            end
        end
    end
end