function copy_strf_files_dir
%copy strf072s cortex files from server to local folders
animal='jack';
area='v1_2';
sessions = main_raw_sessions(animal,area);

for i=1:length(sessions)
    session=sessions(i);
    extNum=1;
    strfFileName=['strf072s.',num2str(extNum)];
    strfFileSource=['V:\thielelab\Groups\ThieleGroup\monkey_data\Jack\_jackgrid\',num2str(session),'\',strfFileName];
    strfFileDestination=['F:\jack\',area,'_strf_files\',num2str(session)];
    if ~exist(strfFileDestination,'dir')
        mkdir(strfFileDestination);
    end
    strfFileDestinationPath=[strfFileDestination,'\strf072s.1'];
    copyfile(strfFileSource,strfFileDestinationPath);
end