function bj_SE_roving_sample_ROC(animal,area,ch,session,sampleContrast,testContrasts,epoch)
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
minusSpontan=0;
if minusSpontan==1
    subfolder=['roving_sample_ROC',num2str(epoch),'_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
elseif minusSpontan==0
    subfolder=['roving_sample_ROC',num2str(epoch),'_wspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
end
sample_act=[];
for cond=1:1:length(testContrasts)
    sample=[];
    %calculate sample-induced activity:
    if epoch==2;
        for n=1:length(matarray{cond,epoch})
            spikeTimes=find(matarray{cond,epoch}{n}<=epochTimes(epoch+1));
            spikeTimes=find(matarray{cond,epoch}{n}(spikeTimes)>epochTimes(epoch));
            sample(n,1)=length(spikeTimes);%tally spikes across trials
        end
        sample_act=[sample_act;sample];%combine across test contrast conditions
    end
end

matSampleName=[num2str(ch),'_',num2str(session),'_',num2str(sampleContrast),'_',area,'_sample_vals'];
matSampleFolder=fullfile('F:','PL','roving_sample_ROC',animal,subfolder);
if ~exist(matSampleFolder,'dir')
    mkdir(matSampleFolder)
end
matSamplePath=fullfile(matSampleFolder,matSampleName);
saveText=['save ',matSamplePath,' sample_act'];
eval(saveText);
