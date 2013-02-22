function fname = spikeData_finder(animal, session, channel,sampleText)

animal = lower(animal);

if ispc
%     if strncmp(getComputerName(),'alex40',6)
        root_roc_vals_folder = fullfile('F:','PL','spikeData',animal);
%     end
end

fname = fullfile( ...
    root_roc_vals_folder , ...
    [num2str(channel),'_',num2str(session),'_',num2str(sampleText),'.mat']);
end