function copy_SessionStandards_files(animal,area)
%Looks up list of representative sessions from which spontaneous activity
%levels will be calculated and used as standards for all other sessions
%for a given channel. Copies representative .nse file from directory
%containing all .nse files across sessions for a given channel, to another
%directory (e.g. in an external HD).

animal = lower(animal);
area = lower(area);

if strcmp(animal,'blanco')
    if strncmp(area,'v4',2)
        allChannels         = [  1   2   3   4   7  12  13  14  15  18  20  22  24  33  34  36  37  38  40  42  49  50  51  52  53  54  55  57  59  60];
        allSessionStandards = [333 333 333 331 333 333 331 330 330 331 331 332 333 335 334 331 330 333 330 332 332 333 332 333 330 334 335 332 331 334];
        %         testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
        newfolder='N:\blanco\v4_sorted_spikes'
        oldfolder='F:\blanco\before_300113\v4_1_sorted_spikes'
    elseif strncmp(area,'v1',2)
        allChannels         = [  8   9  10  11  15  17  19  21  23  25  26  27  28  29  31  44  45  46  48  61  62  63  64];
        allSessionStandards = [351 349 350 350 351 350 349 350 349 350 349 350 350 349 350 351 351 350 349 350 349 350 350];
        %         testContrast=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
        newfolder='N:\blanco\v1_sorted_spikes'
        oldfolder='F:\blanco\before_300113\v1_sorted_spikes'
    end
elseif strcmp(animal,'jack')
    if strncmp(area,'v4',2)
        allChannels         = [ 1  2  3  4  5  6  8 10 24 35 37 39 40 41 49 50 52 53 54 56];
        allSessionStandards = [40 41 41 41 42 40 40 40 41 39 40 40 40 39 39 41 42 41 41 42];
        %         testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
        newfolder='N:\jack\v4_sorted_spikes'
        oldfolder='F:\jack\before_300113\j_v4_1_sorted_spikes'
    elseif strncmp(area,'v1',2)
        allChannels         = [ 7  9 11 12 13 14 15 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32 51 55];
        allSessionStandards = [63 63 64 64 62 64 62 62 63 62 63 63 62 63 64 64 62 62 64 62 64 62 62 62 63];
        %         testContrast=[5 10 15 20 22 25 28 32 35 40 45 50 60 90];
        newfolder='N:\jack\v1_sorted_spikes'
        oldfolder='F:\jack\before_300113\j_v1_sorted_spikes'
    end
end

for i=1:length(allChannels)
    newfolderName=fullfile(newfolder,num2str(allChannels(i)),num2str(allSessionStandards(i)));
    oldfolderName=fullfile(oldfolder,num2str(allChannels(i)),num2str(allSessionStandards(i)));
    mkdir(newfolderName);
    listing=dir(oldfolderName);
    if size(listing,1)>2
        for q=3:size(listing,1)
            fileName=listing(q,1).name;
            oldfileName=fullfile(oldfolderName,fileName);
            newfileName=fullfile(newfolderName,fileName);
            copyfile(oldfileName,newfileName);
        end
    end
end