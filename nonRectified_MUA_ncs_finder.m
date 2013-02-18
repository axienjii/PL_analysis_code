%SPK_NCS_FINDER  Finds Spike NCS files in a known filing system
%   FNAME = NCS_FINDER(ANIMAL, SESSION, CHANNEL)
%   At present, this will work on Newcastle's high-speed PC and Scott's PC.
%   Needs to be extended to work for Alex's PC or Edinburgh cluster, Eddie.

function fname = nonRectified_MUA_ncs_finder(animal, session, channel)


manual_fallback = 0;


animal = lower(animal);

NCSname = ['nonRectified_MUA_Ch_',num2str(channel),'.ncs'];

% try
%     sessfoldername = ncs_folder_finder(animal, session);
% catch ME
%     if strcmp(ME.identifier,'NCSFolderFinder:NoFolder')
%         % Bail now if folder not present
%         % Ask user to locate manually
%         disp('NCS folder not found');
%         fname = manual_locate_ncs(animal, session, channel);
%         return;
%     else
%         rethrow(ME);
%     end 
% end

try
    sessfoldername = nonRectified_MUA_ncs_folder_finder(animal, session, 'spk');
catch ME
    if isunix && ~exist('/media/PLETHRON/','dir') && getusemydrive()
        fname = find_on_drive(animal, session, channel);
        return;
    else
        rethrow(ME);
    end
end
if ispc
        fname = fullfile(sessfoldername,NCSname);
%     if strncmp(getComputerName(),'ALEX44',6)||strncmp(getComputerName(),'ALEX40',6)
%         fname = fullfile(sessfoldername,NCSname);
%     else
%         fname = fullfile(sessfoldername,NCSname);
%     end
else
    fname = fullfile(sessfoldername,'SPK_NCS',NCSname);
end

% Bail now if folder not present
% Ask user to locate manually
if manual_fallback && ~exist(fname,'file')
    fprintf('File %s not found\n',fname);
    fname = manual_locate_ncs(animal, session, channel);
    return;
end

end

% Asks the user to locate the file. Only accepts appropriate file types.
function fname = manual_locate_ncs(animal, session, channel)
    
    instruction_str = ['Locate the NCS file for '...
        capitalise(animal) ' channel ' num2str(channel) ', session ' num2str(session)];
    
    [filename, pathname] = uigetfile('*.ncs;*.Ncs;*.NCS', instruction_str);
    if isequal(filename,0)
       disp('User selected Cancel')
       error('Couldn''t find the NCS file.');
    else
       fname = fullfile(pathname, filename);
       disp(['User selected ', fname])
    end

end

% Find on mydrive, SPARROWHAWK
function fname = find_on_drive(animal, session, channel)
    if isequal({animal, session, channel}, {'jack',72,12})
        fname = '/media/SPARROWHAWK/Documents/summer_project/wave_clus_analysis/jack_V1/session_72/spk_CSC12.ncs';
        if exist(fname,'file'); return; end
    elseif isequal({animal, session, channel}, {'jack',72,18})
        fname = '/media/SPARROWHAWK/Documents/summer_project/wave_clus_analysis/jack_V1/session_72/spk_CSC18.ncs';
        if exist(fname,'file'); return; end
    elseif isequal({animal, session, channel}, {'jack',49,37})
        fname = '/media/SPARROWHAWK/Documents/summer_project/wave_clus_analysis/jack_V4/session_49/spk_CSC37.ncs';
        if exist(fname,'file'); return; end
    elseif isequal({animal, channel}, {'blanco',4})
        fname = ['/media/SPARROWHAWK/Documents/summer_project/raw_files/ncs_files/blanco_v4/4/' num2str(session) '/spk_CSC4.ncs'];
        if exist(fname,'file'); return; end
    end
    error('File does not exist on this drive: /media/SPARROWHAWK.');
end