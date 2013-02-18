function sessions = sessions_available_par(animal, area, expt_type, channel)
% Lists all the sessions with raw nse data or cached testtimes available
% for a given animal and area. If the channel is unspecified, finds the
% raw sessions for the first channel present and assumes they are all the same.

    % default values
    if nargin<3
        expt_type=1;
    end
    if nargin<4
        channel=NaN;
    end
    
    if isnan(channel)
        sessions = raw_sessions_available_par(animal, area, expt_type, channel);
        return;
    end
    
    % Try in case raw files unavailable
    try
        raw_sess = raw_sessions_available_par(animal, area, expt_type, channel);
    catch ME
        disp(ME.message);
        raw_sess = [];
    end
    cached_sess = cached_sessions_available(animal, area, channel);
    sessions = union(raw_sess, cached_sess);
    
end

function sessions = raw_sessions_available_par(animal, area, expt_type, channel)
% Lists all the sessions available for a given animal and area. If the
% channel is unspecified, finds the sessions for the first channel present
% and assumes they are all the same.
    
    % Get the containing folder
    folder = nse_folderfinder_par(animal, area, NaN, channel, expt_type);
    
    % Find a channel if none specified
    if isnan(channel)
        dircontents = dir(folder);
        if length(dircontents)<3
            error('Error with folder: %s\n',folder);
        end
        for ifile=1:length(dircontents)
            file = dircontents(ifile);
            % Ignore non-directories
            if ~file.isdir; continue; end
            % Ignore . and .. and anything else hidden
            fname = file.name;
            if strcmp(fname(1),'.'); continue; end
            ch_fname=fname;
            break;
        end
        folder = fullfile(folder,ch_fname);
    end
    
    % Look through folder for suitable subfolders
    sessions = check_subfolders(folder);
    
    sessions = sort(sessions);
end

function sessions = cached_sessions_available(animal, area, channel)
% Lists all the cached sessions available for a given animal, area and
% channel
    
    % Look up the directory structure
    folder = cachedirname('sptimes_flat', animal, area, channel);
    
    % Look through folder for suitable subfolders
    sessions = check_subfolders(folder);
    
    sessions = sort(sessions);
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
        % Convert underscore to a decimal point, then swap string into num
        session_num = str2double(regexprep(fname,'_','\.'));
        % If we made a number and not NaN, it's a valid sessions number
        if ~isnan(session_num)
            sessions = [sessions session_num];
        end
    end
end
