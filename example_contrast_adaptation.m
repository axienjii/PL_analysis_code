function example_contrast_adaptation
%Plot example figures for channels where noticeable sensory adaptation
%occurred for test-evoked response.

channels=[10 52 53];
sessions=[77 75 46];
% channel=53;
% session=46;
% channel=10;
% session=77;
% figSess=figure('Color',[1,1,1],'Units','Normalized','Position',[0.05, 0.1, 0.3, 0.3]); %
% set(figSess, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
sampleContrast=30;
animal='jack';
figSess=figure('Color',[1,1,1],'Units','Normalized','Position',[0.05, 0.1, 0.3, 0.8]); %
set(figSess, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
for chInd=1:length(channels)
    channel=channels(chInd);
    session=sessions(chInd);
    subplot(3,1,chInd);
    spikeDataName=[num2str(channel),'_',num2str(floor(session)),'_',num2str(sampleContrast)];
    spikeDataFolder=fullfile('F:','PL','spikeData',animal,spikeDataName);
    loadText=['load ',spikeDataFolder,'.mat matarray'];
    eval(loadText);
    
    numBins=100;
    binWidth=512/numBins;
    bins1=0:binWidth:512;
    sigma=3;
    mu=0;
    condCol(8:10,:)=[1 0 0;0.05 0.73 0.24;0 0 1];
    testContrast(8:10)=[31 32 33];
    lineCond=[{'-'} {':'}];
    totalTrials=size(matarray{8,2},1)+size(matarray{9,2},1)+size(matarray{10,2},1);
    histo=zeros(totalTrials,length(bins1));
    for cond=8:10%sample activity, combined across 3 conditions
        for n=1:size(matarray{cond,2},1)
            spikeTimes=matarray{cond,2}{n};%sample
            [N X]=histc(spikeTimes,bins1);%2nd input arg is bin edges
            if size(N,2)==1
                N=N';
            end
            histo(n,1:length(N))=histo(n,1:length(N))+N;
        end
    end
    histostd=std(histo*1000/binWidth,0,1)/totalTrials;
    % histostd=std(histo*1000/binWidth,0,1);
    histo=sum(histo,1)/totalTrials*1000/binWidth;
    fit_val=(-3*sigma):1:(3*sigma);
    y=ones(1,length(fit_val));
    y=(y*(1/(sigma*sqrt(2*pi))).*exp(-(((fit_val-mu).^2)/(2*sigma*sigma))));
    fitted_vector=filter2(y,histo);
    midBins=bins1(2:end)-512/(numBins*2);
    plot(midBins,fitted_vector(1:end-1),'k-','LineWidth',2);hold on%last bin contains number of spikes that equal to the final bin edge- can exclude it
    plot(midBins,fitted_vector(1:end-1)+histostd(1:end-1),'k--');hold on%last bin contains number of spikes that equal to the final bin edge- can exclude it
    plot(midBins,fitted_vector(1:end-1)-histostd(1:end-1),'k--');hold on%last bin contains number of spikes that equal to the final bin edge- can exclude it
    
    bins=512*2:binWidth:512*3;%test act, for each of 3 conditions
    for cond=8:10
        histo=zeros(size(matarray{cond,4},1),length(bins));
        for n=1:size(matarray{cond,4},1)
            spikeTimes=matarray{cond,4}{n};%sample
            [N X]=histc(spikeTimes,bins);%2nd input arg is bin edges
            if size(N,2)==1
                N=N';
            end
            histo(n,1:length(N))=histo(n,1:length(N))+N;
        end
        histostd=std(histo*1000/binWidth,0,1)/size(matarray{cond,4},1);
        %     histostd=std(histo*1000/binWidth,0,1);
        histo=sum(histo,1)/size(matarray{cond,4},1)*1000/binWidth;
        fit_val=(-3*sigma):1:(3*sigma);
        y=ones(1,length(fit_val));
        y=(y*(1/(sigma*sqrt(2*pi))).*exp(-(((fit_val-mu).^2)/(2*sigma*sigma))));
        fitted_vector=filter2(y,histo);
        midBins=bins1(2:end)-512/(numBins*2);
        plot(midBins,fitted_vector(1:end-1),'Color',condCol(cond,:),'LineWidth',2);hold on%last bin contains number of spikes that equal to the final bin edge- can exclude it
        hold on
        plot(midBins,fitted_vector(1:end-1)+histostd(1:end-1),'Color',condCol(cond,:),'LineStyle','--');hold on%last bin contains number of spikes that equal to the final bin edge- can exclude it
        plot(midBins,fitted_vector(1:end-1)-histostd(1:end-1),'Color',condCol(cond,:),'LineStyle','--');hold on%last bin contains number of spikes that equal to the final bin edge- can exclude it
        xlim([midBins(1) midBins(end)]);
        set(gca,'XTick',[midBins(1) midBins(end)],'XTickLabel',[0 512]);
        if chInd==3
            ptext=sprintf('- %d%% test',testContrast(cond));
            text('Position',[360 2.7*cond],'FontSize',14,'String',ptext,'Color',condCol(cond,:));
        end
    end
    if chInd==3
        ptext=sprintf('- 30%% sample',testContrast(cond));
        text('Position',[360 2.7*7],'FontSize',14,'String',ptext,'Color','k');
    end
end
