function bj_plot_example_distinctive_PSTHs
animal='blanco';
area='v4_1';
sessionNums = main_raw_sessions_final(animal,area,[],0);
exampleChannels=[40 7 12 18 49 57];
minusSpontan=0;
sigma=8;
sigma2=1;
sampleContrast=30;
if minusSpontan==1
    subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
elseif minusSpontan==0
    subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_wspontan_',area];%folder for stimulus-evoked responses without any subtraction of spontaneous activity levels
end
figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.5, 0, 0.2, 0.9]);
for i=1:length(exampleChannels)
    meanPSTHact=[];
    bootPSTHact=[];
    stdPSTHact=[];
    %consolidate the PSTHs across sessions
    for j=1:length(sessionNums)
        matPSTHName=[num2str(exampleChannels(i)),'_',num2str(sessionNums(j)),'_',num2str(sampleContrast),'_',area,'_PSTHact'];
        matPSTHPath=fullfile('F:','PL','xcorr',animal,subfolder,matPSTHName);
        loadText=['load ',matPSTHPath,' PSTHact'];
        eval(loadText);
        allPSTHact(j,:)=PSTHact;
    end
    %create a grand mean PSTH
    meanPSTHact(i,:)=mean(allPSTHact,1);
    bootPSTHact(i,:)=gaussfit(sigma,0,meanPSTHact(i,:));%smoothing done for each bootstrapped PSTH
    stdPSTHact(i,:)=std(allPSTHact,0,1);
    subplot(length(exampleChannels),1,i);
    plot(1:length(bootPSTHact(i,:)),meanPSTHact(i,:),'Color','k','LineWidth',2);hold on
    plot(1:length(bootPSTHact(i,:)),meanPSTHact(i,:)-stdPSTHact(i,:),'k:','LineWidth',1);hold on
    plot(1:length(bootPSTHact(i,:)),meanPSTHact(i,:)+stdPSTHact(i,:),'k:','LineWidth',1);hold on
    [yMax timeInd]=max(bootPSTHact(i,:)+stdPSTHact(i,:));
    ylim([0 yMax]);
    xlim([0 512+400]);
%     set(gca,'YTick',[0 yMax]);
%     set(gca,'YTickLabel',[0 round(yMax)]);
    line([timeInd timeInd],[0 yMax],'Color','r','LineStyle',':','LineWidth',1.5);
    title(num2str(exampleChannels(i)));
end
subplot(6,1,6);
xlabel('Time from test stimulus onset (ms)');
ylabel('Firing rate (spikes/s)');
imageName=['example_channels_distinctive_waveforms_',animal,'_',area];
imageFolder=fullfile('F:','PL','xcorr',animal,subfolder);
imagePath=fullfile(imageFolder,imageName);
printtext=['print -dpng -r300 ',imagePath];
set(gcf,'PaperPositionMode','auto')
eval(printtext);
close all