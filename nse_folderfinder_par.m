function folder = nse_folderfinder_par(animal,area,session,channel,expt_type)
% Function to navigate file structure for finding NSE files
    
    % default values
    if nargin<4
        channel=NaN;
    end
    if nargin<3
        session=NaN;
    end
    if ~isnan(session)
        [ignore,ignore,ignore,expt_type] = session_metadata(session, animal);
    elseif nargin<5
        expt_type = 1;
    end
    
    animal = lower(animal);
    area = lower(area);
    
    % OS dependent addresses
    if ispc && ~getusemydrive()
%         if strncmp(getComputerName(),'ALEX44',6)
%             root_dir = 'C:';
%         elseif strncmp(getComputerName(),'ALEX40',6)
            root_dir = 'F:';
%         end
    else
        root_dir = fullfile(getrootdir(),'raw_files');
    end
    
    if strcmp(area,'v4')
        regionphase = [area,'_',num2str(expt_type)];
    elseif strcmp(area,'v1') && expt_type~=1
        regionphase = [area,'_',num2str(expt_type)];
    else
        regionphase = area;
    end
    
    if strcmp(animal,'jack')
        subdirname = ['j_',regionphase,'_sorted_spikes'];
    else
        subdirname = [regionphase,'_sorted_spikes'];
    end
    
    folder = fullfile(root_dir,animal,subdirname);
    
    if ~isnan(channel)
        folder = fullfile(folder,num2str(channel));
        if ~isnan(session)
            % Need to swap decimals back to underscores
            folder = fullfile(folder,regexprep(num2str(session),'\.','_'));
        end
    end
end
