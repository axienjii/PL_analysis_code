function move_images_to_folder(animal,area)
%Moves saved image files into folder, categorised by file format (png or
%eps).

animal = lower(animal);
area = lower(area);

if strcmp(animal,'blanco')
    if strncmp(area,'v4',2)
        allChannels         = [  1   2   3   4   7  12  13  14  15  18  20  22  24  33  34  36  37  38  40  42  49  50  51  52  53  54  55  57  59  60];
    elseif strncmp(area,'v1',2)
        allChannels         = [  8   9  10  11  15  17  19  21  23  25  26  27  28  29  31  44  45  46  48  61  62  63  64];
    end
elseif strcmp(animal,'jack')
    if strncmp(area,'v4',2)
        allChannels         = [ 1  2  3  4  5  6  8 10 24 35 37 39 40 41 49 50 52 53 54 56];
    elseif strncmp(area,'v1',2)
        allChannels         = [ 7  9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
    end
end

sessions = main_raw_sessions(animal,area);

formats=[{'eps'} {'png'}];
for i=1:length(allChannels)
    for j=1:2
        newfolderName=fullfile('F:','PL','PSTHs',animal,num2str(allChannels(i)),formats{j});
        mkdir(newfolderName);
        for k=1:length(sessions)
            fileName=[num2str(allChannels(i)),'_',num2str(sessions(k)),'_30.',formats{j}];
            oldPathName=fullfile('F:','PL','PSTHs',animal,num2str(allChannels(i)),fileName);
            newPathName=fullfile(newfolderName,fileName);
            if exist(oldPathName,'file')
                movefile(oldPathName,newPathName);
            end
        end
    end
end