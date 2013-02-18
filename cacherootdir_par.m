function folder = cacherootdir()

if strncmp(getComputerName(),'eddie',5)
    % Cache and raw files stored seperately on Eddie
    root_dir = '/exports/work/inf_ndtc/s1145806';
else
    root_dir = getrootdir_par();
end
folder = fullfile(      ...
    root_dir            , ...
    'generated_files'   , ...
    'cache'             );
end