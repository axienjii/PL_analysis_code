function move_artifact_images
%written on 20/02/13 by Xing to automate moving of files from artifact-free folders (to which
%they were accidentally written at first) to artifact-present folders.
animal='jack';
area='v1_2';
fileFormats=[{'eps'} {'png'}];
sampleContrasts=[20 30 40];
channels = main_channels(animal,area);
sessions = main_raw_sessions(animal,area);
for j=1:length(channels)
    channel=channels(j);
    for i=1:length(sessions)
        session=sessions(i);
        for k=1:2
            for m=1:3
                sampleContrast=sampleContrasts(m);
                sourcePath=['F:\PL\PSTHs\jack\',num2str(channel),'\',fileFormats{k},'\',num2str(channel),'_',num2str(session),'_',num2str(sampleContrast),'.',fileFormats{k}];
                destinPath=['F:\PL\with_artifact_trials\PSTHs\jack\',num2str(channel),'\',fileFormats{k},'\',num2str(channel),'_',num2str(session),'_',num2str(sampleContrast),'.',fileFormats{k}];
                movefile(sourcePath,destinPath);
            end
        end
    end
end
    