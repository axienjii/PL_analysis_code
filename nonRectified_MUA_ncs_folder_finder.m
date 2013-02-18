%NCS_FOLDER_FINDER Finds the root NCS folder on this system
function foldername = nonRectified_MUA_ncs_folder_finder(animal, session, type)

if nargin<3
    type = 'spk';
end

% Lookup hostname of computer. Different computers have files in different
% locations.
% Can switch based on this
% compname = getComputerName();

% This list needs to be extended for other pcs
if ispc
        root_ncs_folder = fullfile('N:', animal);
%         root_ncs_folder = fullfile('I:', animal, 'MehdiCheetahData');
%     if strncmp(getComputerName(),'ALEX44',6)
%         if strcmpi(animal,'blanco')
%             root_ncs_folder = fullfile('C:', animal);
%         elseif strcmpi(animal,'jack')
%             root_ncs_folder = fullfile('E:', animal);
%         end
% %         root_ncs_folder = fullfile('C:', animal, 'OriginalRange');
%     elseif strncmp(getComputerName(),'ALEX40',6)
%         root_ncs_folder = fullfile('I:', animal, 'MehdiCheetahData');
%     else
%         root_ncs_folder = fullfile('I:', animal, 'MehdiCheetahData');
%     end
else
    if strncmp(getComputerName(),'eddie',5)
        root_ncs_folder = fullfile('/exports/work/scratch/s1145806/pl_semiraw_SPK', animal, 'MehdiCheetahData');
    elseif exist('/media/PLETHRON/','dir')
        if strcmpi(type,'spk')
            root_ncs_folder = fullfile('/media/PLETHRON/','pl_semiraw_SPK', animal, 'MehdiCheetahData');
        else
            root_ncs_folder = fullfile('/media/PLETHRON/','pl_semiraw', animal, 'MehdiCheetahData');
        end
    elseif exist('/media/SPARROWHAWK/','dir')
        root_ncs_folder = fullfile('/media/SPARROWHAWK/Documents/summer_project/raw_files/ncs_files', animal, 'MehdiCheetahData');
    else
        error('Location unknown');
    end
end

% Bail now if folder not present
if ~exist(root_ncs_folder,'dir')
    ME = MException('NCSFolderFinder:NoRootFolder', ...
       'No NCS folder found in: %s',root_ncs_folder);
    throw(ME);
end

% Find the folder for this session
target_sub_folder = num2str(session);
listing=dir(root_ncs_folder);
if size(listing,1)>2
    for ifile = 1:length(listing)
        file = listing(ifile);
        if ~file.isdir
            continue;
        end
        file_name = file.name;
        % Can find either the folder with the same name as the session num
        % or the name followed by an underscore and then some more stuff
        if ( strcmp(target_sub_folder, file_name) || ...
                strncmp([target_sub_folder,'_'], file_name, length(target_sub_folder)+1))
            foldername = fullfile(root_ncs_folder,file_name);
            return;
        end
    end
end
ME = MException('NCSFolderFinder:NoFolder', ...
   'No NCS folder found in: %s',fullfile(root_ncs_folder,target_sub_folder));
throw(ME);
end