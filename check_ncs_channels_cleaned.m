function check_ncs_channels_cleaned(animal,area,sessions)
%Checks whether artifact-free .ncs files have been generated. Prints list
%of missing channels for each session of interest, to screen.
channels = main_channels(animal,area);
missin=[];present=0;
for j=1:length(sessions)
    missin{j}=sessions(j);
    for i=1:length(channels)
        if strcmp(animal,'blanco')
            NCSfilename=['C:\pl_spk_artifactRemoved\blanco\MehdiCheetahData\',num2str(sessions(j)),'\SPK_NCS_artifactRemoved\spk_CSC',num2str(channels(i)),'_artifactRemoved.ncs'];
        elseif strcmp(animal,'jack')
            NCSfilename=['E:\pl_spk_artifactRemoved\jack\MehdiCheetahData\',num2str(sessions(j)),'\SPK_NCS_artifactRemoved\spk_CSC',num2str(channels(i)),'_artifactRemoved.ncs'];
        end
        if exist(NCSfilename,'file')
            present=1;
        end
        if present==0
            missin{j}=[missin{j} channels(i)];
        end
        present=0;
    end
    missin{j}
end