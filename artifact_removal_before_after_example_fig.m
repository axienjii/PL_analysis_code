function artifact_removal_before_after_example_fig
fighandle1=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);

spikeRowInd=1;
load('F:\PL\old version\spikeData\blanco\4_332_30.mat')
matarrayArt=matarray;
load('F:\PL\spikeData\blanco\plotRedArtifacts\correct_trials_only\4_332_30.mat')
subplot(1,2,1);
for condInd=1:size(matarrayArt,1)
    for trialInd=1:size(matarrayArt{condInd,3})
        for spikeInd=1:length(matarrayArt{condInd,3}{trialInd})
            plot([matarrayArt{condInd,3}{trialInd}(spikeInd) matarrayArt{condInd,3}{trialInd}(spikeInd)],[spikeRowInd-0.5 spikeRowInd+0.5],'Color','k');hold on
            spikeRowInd=spikeRowInd+1;
        end
    end
end
title('before artifact removal');
ylabel('trial number');
xlabel('time (  )');
spikeRowInd=1;
subplot(1,2,2);
for condInd=1:size(matarray,1)
    for trialInd=1:size(matarray{condInd,3})
        for spikeInd=1:length(matarray{condInd,3}{trialInd})
            plot([matarray{condInd,3}{trialInd}(spikeInd) matarray{condInd,3}{trialInd}(spikeInd)],[spikeRowInd-0.5 spikeRowInd+0.5],'Color','k');hold on
            spikeRowInd=spikeRowInd+1;
        end
    end
end
title('after artifact removal');