function root_dir = getrootdir_par()

    [usemydrive drive_letter] = getusemydrive();

    if ispc
        if usemydrive
            root_dir = fullfile([drive_letter ':'],'Documents','summer_project');
        else
%             if strncmpi(getComputerName(),'ALEX44',6)
%                 root_dir = 'C:';
%             elseif strncmpi(getComputerName(),'ALEX40',6)
                root_dir = 'F:';
%             end
        end
    elseif isunix
        if exist('/media/PLETHRON','dir')
            root_dir = fullfile('/media/PLETHRON');
        elseif usemydrive
            root_dir = fullfile('/media/SPARROWHAWK','Documents','summer_project');
        elseif exist('/disk/scratch','dir')
            root_dir = fullfile('/disk/scratch','Documents','summer_project');
        elseif strncmp(getComputerName(),'eddie',5)
            root_dir = '/exports/work/scratch/s1145806';
        else
            root_dir = fullfile(getenv('HOME'),'Documents','summer_project');
        end
    end
end