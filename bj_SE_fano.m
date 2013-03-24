function bj_SE_fano(animal,area,ch,session,sampleContrast,testContrasts,epoch)
%Written by Xing 09/09/10
%Modified from use_time_periods2, calculates PSTH values during test
%presentation at highest contrast, writes activity levels to file, PSTHact. Activity
%compiled across sessions for each cell, for cross correlation analysis
%later, by function read_blanco_SE_xcorr.
%
%Naming of folders: PSTH45_images changed to e.g.
%-PSTH45_images_sm6ms_wspontan: means that epochs 4 & 5 are combined into a
%single array called PSTHact, gaussfit is carried out with sigma of 3 ms
%(to produce smoothing over 6 ms), and spontaneous activity levels are not
%subtracted from stimulus-evoked responses. Smoothing of 6 ms, with spontan.    
%-PSTH45_images_sm10ms_mpontan: means that gaussfit has sigma of 5 ms
%(smoothing over 10 ms), and spontaneous activity has been subtracted from
%PSTHact. Smoothing of 10 ms, minus spontan.  

epochTimes=[-512 0 512 1024 1536 1936]; 
matName=[num2str(ch),'_',num2str(session),'_',num2str(sampleContrast)];
matPath=fullfile('F:','PL','spikeData',animal,matName);
loadText=['load ',matPath,' matarray'];
eval(loadText);
binwidth=1;%1 ms
minusSpontan=0;
if minusSpontan==1
    subfolder=['fano',num2str(epoch),'_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
elseif minusSpontan==0
    subfolder=['fano',num2str(epoch),'_wspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
end
for cond=1:1:length(testContrasts)
    %calculate spontaneous activity:
    if epoch==1;
        bins=epochTimes(epoch)+binwidth/2:binwidth:epochTimes(epoch+1)-binwidth/2;
        spontan1=zeros(1,length(bins));
        for n=1:length(matarray{cond,epoch})
            spikeTimes=find(matarray{cond,epoch}{n}<=epochTimes(epoch+1));
            spikeTimes=find(matarray{cond,epoch}{n}(spikeTimes)>epochTimes(epoch));%time stamps from -150 to 0 ms relative to sample onset
            spontan1(n)=length(spikeTimes);%tally number of spikes across trials
        end
    end
    %calculate PSTH for epochs 4 & 5 combined
    if epoch==4;
        bins4=epochTimes(epoch)+binwidth/2:binwidth:epochTimes(epoch)+512-binwidth/2;% to 512 ms
        for n=1:length(matarray{cond,epoch})
            spikeTimes=find(matarray{cond,epoch}{n}<=epochTimes(epoch+1));
            spikeTimes=find(matarray{cond,epoch}{n}(spikeTimes)>epochTimes(epoch));
            test(n)=length(spikeTimes);%tally spikes across trials
        end
        bins=bins4;%just calculate fano factor based on epoch 4
        fano_vals(cond)=var(test)/mean(test);
    end
    if epoch==5;
        bins5=epochTimes(epoch)+binwidth/2:binwidth:epochTimes(epoch)+400-binwidth/2;% to 512 ms
        postTest=zeros(1,length(bins5));
        for n=1:length(matarray{cond,epoch})
            spikeTimes=find(matarray{cond,epoch}{n}<=epochTimes(epoch+1));
            spikeTimes=find(matarray{cond,epoch}{n}(spikeTimes)>epochTimes(epoch));
            postTest(n)=length(spikeTimes);%tally spikes across trials
        end
    end
    %     bins=[bins4 bins5];
end

matFanoName=[num2str(ch),'_',num2str(session),'_',num2str(sampleContrast),'_',area,'_fano_vals'];
matFanoFolder=fullfile('F:','PL','fano',animal,subfolder);
if ~exist(matFanoFolder,'dir')
    mkdir(matFanoFolder)
end
matFanoPath=fullfile(matFanoFolder,matFanoName);
saveText=['save ',matFanoPath,' fano_vals'];
eval(saveText);
