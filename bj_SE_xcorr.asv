function bj_SE_xcorr(animal,area,ch,session,sigma,sampleContrast,testContrasts)
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

minusSpontan=0;
epochTimes=[-512 0 512 1024 1536 1936]; 
cond=length(testContrasts);%either 14 or 12- calculate cross-correlation values based on response elicited by highest contrast stimulus
matName=[num2str(ch),'_',num2str(session),'_',num2str(sampleContrast)];
matPath=fullfile('F:','PL','spikeData',animal,matName);
loadText=['load ',matPath,' matarray'];
eval(loadText);
binwidth=1;%1 ms
%calculate spontaneous activity:
epoch=1;
bins=epochTimes(epoch)+binwidth/2:binwidth:epochTimes(epoch+1)-binwidth/2;
spontan1=zeros(1,length(bins));
for n=1:length(matarray{cond,epoch})
    spikeTimes=find(matarray{cond,epoch}{n}<=epochTimes(epoch+1));
    spikeTimes=find(matarray{cond,epoch}{n}(spikeTimes)>epochTimes(epoch));%time stamps from -150 to 0 ms relative to sample onset
    [N X]=hist(matarray{cond,epoch}{n}(spikeTimes),bins);
    if size(N,2)==1
        N=N';
    end
    spontan1(1:length(N))=spontan1(1:length(N))+N;%tally spikes across trials
end
spontan1(1,:)=spontan1(1,:)*1000/(binwidth*length(matarray{cond,epoch}));%average activity per ms
spontan1(1,:)=gaussfit(sigma,0,spontan1(1,:));
mean_spontan1=mean(spontan1);
std_spontan1=std(spontan1);

%calculate PSTH for epochs 4 & 5 combined
epoch=4;
bins4=epochTimes(epoch)+binwidth/2:binwidth:epochTimes(epoch)+512-binwidth/2;% to 512 ms
test=zeros(1,length(bins4));
for n=1:length(matarray{cond,epoch})
    spikeTimes=find(matarray{cond,epoch}{n}<=epochTimes(epoch+1));
    spikeTimes=find(matarray{cond,epoch}{n}(spikeTimes)>epochTimes(epoch));
    [N X]=hist(matarray{cond,epoch}{n}(spikeTimes),bins4);
    if size(N,2)==1
        N=N';
    end
    test(1:length(N))=test(1:length(N))+N;%tally spikes across trials
end
test(1,:)=test(1,:)*1000/(binwidth*length(matarray{cond,epoch}));%average activity per ms

epoch=5;
bins5=epochTimes(epoch)+binwidth/2:binwidth:epochTimes(epoch)+400-binwidth/2;% to 512 ms
postTest=zeros(1,length(bins5));
for n=1:length(matarray{cond,epoch})
    spikeTimes=find(matarray{cond,epoch}{n}<=epochTimes(epoch+1));
    spikeTimes=find(matarray{cond,epoch}{n}(spikeTimes)>epochTimes(epoch));
    [N X]=hist(matarray{cond,epoch}{n}(spikeTimes),bins5);
    if size(N,2)==1
        N=N';
    end
    postTest(1:length(N))=postTest(1:length(N))+N;%tally spikes across trials
end
postTest(1,:)=postTest(1,:)*1000/(binwidth*length(matarray{cond,epoch}));%average activity per ms

bins=[bins4 bins5];
PSTHact=[test(1,:) postTest(1,:)];
PSTHact=gaussfit(sigma,0,PSTHact);%smoothing done here
figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.1, 0.3, 0.7, 0.4]);
%to subtract spontan act:
if mean_spontan1<=mean(PSTHact)&&minusSpontan==1
    oriPSTHact=PSTHact;%make copy of original act, without subtraction of spontaneous levels
    for actInd=1:length(PSTHact)
        if PSTHact(actInd)-mean_spontan1>=0
            PSTHact(actInd)=PSTHact(actInd)-mean_spontan1;
        else
            PSTHact(actInd)=0;
        end
    end
    plot(bins,oriPSTHact,'Color','c');hold on
    maxY=[max(oriPSTHact) max(PSTHact)];
    maxY=max(maxY);
    ptext=sprintf('ch %s   session %d     spontan %.1f spikes/s',num2str(ch),session,mean_spontan1);
else
    maxY=max(PSTHact);
    ptext=sprintf('ch %s   session %d',num2str(ch),session);
end
plot(bins,PSTHact,'Color','k')
line([1536 1536],[0 maxY],'LineStyle',':','Color','k');
set(gca,'XLim',[epochTimes(4) epochTimes(4)+512+400]);
set(gca,'XTick',[epochTimes(4) epochTimes(4)+512 epochTimes(4)+512+400]);
set(gca,'XTickLabel',[1024 1536 1936]);
set(gca,'YLim',[0 maxY]);
text('Position',[1024 1.05*maxY],'FontSize',9,'String',ptext);

imageName=[num2str(ch),'_',num2str(session),'_',num2str(sampleContrast),'_',area];
if minusSpontan==1
    subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
elseif minusSpontan==0
    subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_wspontan_',area];%folder for stimulus-evoked responses without any subtraction of spontaneous activity levels
end
imageFolder=fullfile('F:','PL','xcorr',animal,subfolder);
if ~exist(imageFolder,'dir')
    mkdir(imageFolder);
end
imagePath=fullfile(imageFolder,imageName);
printtext=['print -dpng ',imagePath];
set(gcf,'PaperPositionMode','auto')
eval(printtext);

matPSTHName=[num2str(ch),'_',num2str(session),'_',num2str(sampleContrast),'_',area,'_PSTHact'];
matPSTHPath=fullfile('F:','PL','xcorr',animal,subfolder,matPSTHName);
saveText=['save ',matPSTHPath,' PSTHact'];
eval(saveText);

close all
