function fname = nev_finder_par(animal,area,session,varargin)
% Function to navigate file structure for finding NEV files
    
    folder = nev_folderfinder_par(animal,area,session,varargin{:});
    
    animal = lower(animal);
    
    % OS dependent addresses
    if ispc
        if strcmp(animal,'blanco')
            file_name = 'Events.nev';
        elseif strcmp(animal,'jack')
            file_name = 'Events.Nev';
        end
    elseif isunix
        file_name = 'Events.nev';
    end
    
    fname = fullfile(folder,file_name);
    
    if ~exist(fname,'file')
%         fprintf('!!!!!!!!!!! File not found !!!!!!!!!!!!\nNEV for %s %s, session %s\n%s\n',...
%             animal,area,num2str(session),fname);
%         fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        ME = MException('NEVFinder:NoFile', ...
           'NEV file not present: %s',fname);
        throw(ME);
    end
end

function folder = nev_folderfinder_par(animal,area,session,varargin)
% Function to navigate file structure for finding NEV files
    
    % default values
%     if nargin<4
        [ignore, ignore, ignore, expt_type] = session_metadata(session, animal);
%     end
    
    animal = lower(animal);
    area = lower(area);
    
    % OS dependent addresses
    if ispc && ~getusemydrive()
        if strcmp(animal,'blanco')
            if strncmp(getComputerName(),'ALEX44',6)
%                 root_dir = 'C:';
            elseif strncmp(getComputerName(),'alex40',6)
                root_dir = 'F:';
            end
        elseif strcmp(animal,'jack')
            root_dir = 'V:\thielelab\Groups\ThieleGroup';
        end
    else
        root_dir = fullfile(getrootdir(),'raw_files');
        if ~exist(root_dir,'dir')
            root_dir = '/media/PLETHRON/raw_files';
        end
    end
    
    % Different experiment parts
    if strcmp(area,'v4')
        regionphase = [area,'_',num2str(expt_type)];
    elseif strcmp(area,'v1') && expt_type~=1
        regionphase = [area,'_',num2str(expt_type)];
    else
        regionphase = area;
    end
    
    % animal dependent addresses
    if strcmp(animal,'jack')
        capital_animal = animal;
        capital_animal(1)=upper(capital_animal(1));
        folder = fullfile(root_dir,...
            'monkey_data',...
            capital_animal,...
            ['_',animal,'grid'],...
            ['j_',regionphase,'_events_files']...
            );
    elseif strcmp(animal,'blanco')
        folder = fullfile(root_dir,...
            animal,...
            [regionphase,'_events_files']...
            );
    end
    
    % Need to swap decimal point back to underscore
    folder = fullfile(folder,regexprep(num2str(session),'\.','_'));
end