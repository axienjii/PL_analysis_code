function read_ROC_old_new(cutoff,excludeSuppressed,normalised)
analysisType='ROC';
analysisType='ROC_zero_one';
%Written by Xing 19/05/13
%Compare psychometric function parameters between new and old methods of
%ROC/AUROC value generation
onExternalHD=0;
if onExternalHD==1
    rootFolder='G:\PL_backup_060413';
else
    rootFolder='F:';
end
figGauss=figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.45, 0.65]); %
set(figGauss, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 3.305 3.325/0.4*0.45 3.305/5*65]);
animals=[{'blanco'} {'jack'}];
areas=[{'v4_1'} {'v1_1'}];
allTableStats=[];
allParTableStats=[];
mpallTableStats=[]
for areaInd=1:length(areas)
    area=areas{areaInd};
    for animalInd=1:length(animals)
        animal=animals{animalInd};
%         if cutoff==1
%             loadText=['load ',rootFolder,'\PL\',analysisType,'\',animal,'\new_vs_old_sglrocmeanchannels\cumulative_ROCs_old_new_',area,'_30.mat'];
%             loadText=['load ',rootFolder,'\PL\',analysisType,'\',animal,'\new_vs_old_sglroc3acrosschannels\cumulative_ROCs_old_new_',area,'_30.mat'];
%         else
            subFolder='new_vs_old_sglrocmeanchannels';
            subFolder='new_vs_old_sglroc3acrosschannels';
            if excludeSuppressed==1
                subFolder=[subFolder,'_excludeSuppressed'];
            end
            if normalised
                subFolder=[subFolder,'_normalised'];
            end
            loadText=['load ',rootFolder,'\PL\',analysisType,'\',animal,'\',subFolder,'\cumulative_ROCs_old_new_',area,'_30_cutoff',num2str(cutoff*10),'.mat'];  
            %         end
        eval(loadText)
        newValues=[slopeNeuroNew;PNENew;minRateNew;maxRateNew];
        oldValues=[slopeNeuroOld;PNEOld;minRateOld;maxRateOld];
        dfs=[statsS.df statsP.df statsmin.df statsmax.df];
        tvals=[statsS.tstat statsP.tstat statsmin.tstat statsmax.tstat];
        pvals=[pS pP pmin pmax];
        for paramInd=1:4
            subplot(4,4,(paramInd-1)*4+animalInd+(areaInd-1)*2);
            %calculate mean, SD for each parameter, for new & old methods:
            tableStats(paramInd,1)=mean(oldValues(paramInd,:));
            tableStats(paramInd,2)=std(oldValues(paramInd,:));
            tableStats(paramInd,3)=mean(newValues(paramInd,:));
            tableStats(paramInd,4)=std(newValues(paramInd,:));
            %ct-test stats for difference between old & new values:
            tableStats(paramInd,5)=dfs(paramInd);
            tableStats(paramInd,6)=tvals(paramInd);
            tableStats(paramInd,7)=pvals(paramInd);
%             plot(1:length(slopeNeuroNew),newValues(paramInd,:),'ro','LineStyle','none','MarkerFaceColor','r');hold on
%             plot(1:length(slopeNeuroNew),oldValues(paramInd,:),'bo','LineStyle','none','MarkerFaceColor','b');
%             plot(1:length(slopeNeuroNew),newValues(paramInd,:),'rx','LineStyle','none');hold on
%             plot(1:length(slopeNeuroNew),oldValues(paramInd,:),'bx','LineStyle','none','MarkerSize',7);
%             plot(1:length(slopeNeuroNew),oldValues(paramInd,:),'bx','LineStyle','none','MarkerSize',7);hold on
%             plot(1:length(slopeNeuroNew),newValues(paramInd,:),'r+','LineStyle','none');
            plot(1:length(slopeNeuroNew),newValues(paramInd,:),'ro','LineStyle','none');hold on
            plot(1:length(slopeNeuroNew),oldValues(paramInd,:),'b+','LineStyle','none','MarkerSize',7);
%             plot(1:length(slopeNeuroNew),newValues(paramInd,:),'ro','LineStyle','none');hold on
%             plot(1:length(slopeNeuroNew),oldValues(paramInd,:),'bx','LineStyle','none','MarkerSize',7);
            xlim([0 length(slopeNeuroNew)+1]);
            axis square
        end
        %calculate partial correlations for each parameter (based on new data only):
        [rho p]=partialcorr(slopeNeuroNew',(1:length(slopeNeuroNew))',[PNENew' minRateNew' maxRateNew'])
        parTableStats(1,8:10)=[length(slopeNeuroNew)-2 rho p];
        [rho p]=partialcorr(PNENew',(1:length(PNENew))',[slopeNeuroNew' minRateNew' maxRateNew'])
        parTableStats(2,8:10)=[length(PNENew)-2 rho p];
        [rho p]=partialcorr(minRateNew',(1:length(minRateNew))',[slopeNeuroNew' PNENew' maxRateNew'])
        parTableStats(3,8:10)=[length(minRateNew)-2 rho p];
        [rho p]=partialcorr(maxRateNew',(1:length(maxRateNew))',[slopeNeuroNew' PNENew' minRateNew'])
        parTableStats(4,8:10)=[length(maxRateNew)-2 rho p];
        allParTableStats=[allParTableStats;parTableStats];
        %calculate regular correlations
        [rho p]=corr(slopeNeuroNew',(1:length(slopeNeuroNew))')
        tableStats(1,8:10)=[length(slopeNeuroNew)-2 rho p];
        [rho p]=corr(PNENew',(1:length(PNENew))')
        tableStats(2,8:10)=[length(PNENew)-2 rho p];
        [rho p]=corr(minRateNew',(1:length(minRateNew))')
        tableStats(3,8:10)=[length(minRateNew)-2 rho p];
        [rho p]=corr(maxRateNew',(1:length(maxRateNew))')
        tableStats(4,8:10)=[length(maxRateNew)-2 rho p];
        allTableStats=[allTableStats;tableStats];    
        %calculate meaningful partial correlations
        [rho p]=partialcorr(slopeNeuroNew',(1:length(slopeNeuroNew))',[minRateNew' maxRateNew'])
        mptableStats(1,8:10)=[length(slopeNeuroNew)-2 rho p];
        [rho p]=partialcorr(PNENew',(1:length(PNENew))',[minRateNew' maxRateNew'])
        mptableStats(2,8:10)=[length(PNENew)-2 rho p];
        [rho p]=partialcorr(minRateNew',(1:length(minRateNew))',maxRateNew')
        mptableStats(3,8:10)=[length(minRateNew)-2 rho p];
        [rho p]=partialcorr(maxRateNew',(1:length(maxRateNew))',minRateNew')
        mptableStats(4,8:10)=[length(maxRateNew)-2 rho p];
        mpallTableStats=[mpallTableStats;mptableStats];      
    end
end
if normalised==0
    subplot(4,4,1);
    xlabel('contrast (%)');
    ylabel('slope');
    title('monkey 1');
    subplot(4,4,2);
    title('monkey 2');
    subplot(4,4,3);
    title('monkey 1');
    subplot(4,4,4);
    title('monkey 2');
    subplot(4,4,5);
    ylabel('PNE');
    subplot(4,4,9);
    ylabel('minimum');
    subplot(4,4,13);
    ylabel('maximum');
    subplot(4,4,1);
    ylim([0.015 0.05]);
    subplot(4,4,2);
    ylim([0 0.035]);
    subplot(4,4,3);
    ylim([0.007 0.025]);
    subplot(4,4,5);
    ylim([30 39]);
    subplot(4,4,7);
    ylim([28 38]);
    subplot(4,4,12);
    ylim([0 0.4]);
    subplot(4,4,13);
    ylim([0.7 1]);
    subplot(4,4,14);
    ylim([0.7 1]);
    subplot(4,4,15);
    ylim([0.7 1]);
    subplot(4,4,16);
    ylim([0.7 1]);
end
% subplot(4,4,3);
% ylim([0.0025 0.007]);
% subplot(4,4,4);
% ylim([0.011 0.019]);
% subplot(4,4,8);
% ylim([33 40]);
% subplot(4,4,9);
% ylim([-0.4 0.4]);
% subplot(4,4,10);
% ylim([-0.2 0.5]);
% subplot(4,4,11);
% ylim([-0.2 0.4]);
% subplot(4,4,13);
% ylim([0.2 1]);
% subplot(4,4,14);
% ylim([0.2 0.8]);
% subplot(4,4,16);
% ylim([0.78 1]);
% %ignore outliers:
% subplot(4,4,9);
% ylim([0.05 0.4]);
% subplot(4,4,10);
% ylim([0.25 0.45]);
% subplot(4,4,13);
% ylim([0.2 0.6]);
% subplot(4,4,14);
% ylim([0.2 0.5]);
imagename=['bj_old_new_ROC_4_params_sessions_cutoff',num2str(cutoff*10)];
if excludeSuppressed==1
    imagename=['bj_old_new_ROC_4_params_sessions_cutoff',num2str(cutoff*10),'_excludeSuppressed'];
end
pathname=fullfile(rootFolder,'PL','ROC',imagename);
printtext=sprintf('print -dpng %s.png',pathname);
set(gcf,'PaperPositionMode','auto')
eval(printtext);
cutoff
allTableStats