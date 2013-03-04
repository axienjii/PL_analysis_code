function channels = channels_available_par(animal, area, expt_type)
% Lists all the channels/sorted neuron numbers available for analysis for a
% given animal and area
    
    % default values
    if nargin<3
        expt_type=1;
    end
    
    raw_channels = raw_channels_available_par(animal, area, expt_type);
    cached_channels = cached_channels_available(animal, area);
    channels = union(raw_channels, cached_channels);
    
end

function channels = raw_channels_available_par(animal,area,expt_type)
% Lists all the channels/sorted neuron numbers available for analysis for a
% given animal and area
    
    % Get the containing folder
    if ispc
        channels=cached_channels_available(animal,area)
    else
        folder = nse_folderfinder_par(animal,area,NaN,NaN,expt_type);
    
        channels = check_subfolders(folder);
        
        channels = sort(channels);
    end
end

function channels = cached_channels_available(animal,area)
% Lists all the cached sessions available for a given animal, area and
% channel
    
    % Look up the directory structure
    folder = cachedirname('sptimes_flat', animal, area);
    
    % Look through folder for suitable subfolders
    channels = check_subfolders(folder);
    if ispc
        if strcmpi(animal,'blanco')
            if strcmpi(area,'v4')||strcmpi(area,'v4_2')
                channels= [1 2 3 4 7 12 13 14 18 20 22 24 33 34 36 37 38 40 42 49 50 51 52 53 54 55 57 59 60];
            else
                channels=[8 9 10 11 15 17 19 21 23 25 26 27 28 29 31 44 45 46 48 61 62 63 64];
            end
        elseif strcmpi(animal,'jack')
            if strcmpi(area,'v4')||strcmpi(area,'v4_2')||strcmpi(area,'v4_3')
                channels=[1:5 6 8 10 24 35 37 39 40 41 49 50 52:54 56];
            else
                channels=[7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
            end
        end
    end
    channels = sort(channels);
end

function sessions = check_subfolders(folder)
    
    dircontents = dir(folder);
    sessions = [];
    
    % Look through folder for suitable subfolders
    for ifile=1:length(dircontents)
        file = dircontents(ifile);
        % Ignore non-directories
        if ~file.isdir; continue; end
        % Ignore . and .. and anything else hidden
        fname = file.name;
        if strcmp(fname(1),'.'); continue; end
        session_num = str2double(fname);
        if ~isnan(session_num)
            sessions = [sessions session_num];
        end
    end
end
