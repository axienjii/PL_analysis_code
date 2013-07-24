function bj_SE_roving_sample_ROC(animal,area,ch,session,sampleContrast,testContrasts,epoch)
%Written by Xing 09/09/10

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
