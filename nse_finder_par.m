function fname = nse_finder_par(animal,area,session,channel,varargin)
% Function to navigate file structure for finding NSE files
    
    folder = nse_folderfinder_par(animal,area,session,channel,varargin{:});
    
    listing=dir(folder);
    if size(listing,1)>2
        for ifile = 1:length(listing)
            file = listing(ifile);
            if file.isdir
                continue;
            end
            file_name = file.name;
            if strcmpi(file_name(end-3:end),'.nse')
                fname = fullfile(folder,file_name);
                return;
            end
        end
    end
%     fprintf('!!!!!!!!!!!!! File not found !!!!!!!!!!!!!!\nNSE for Channel %s Session %s\n%s\n',...
%         num2str(channel),num2str(session),folder);
%     fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    ME = MException('NSEFinder:NoFile', ...
       'No NSE file found in: %s',folder);
    throw(ME);
end
