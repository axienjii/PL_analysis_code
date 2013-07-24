function bj_earlylate_behav(psychoOnly)
%Modified from plot_behav. Rewritten by Xing 29/10/12 after HD crash
%Generates boxplots with results, comparing first 30% with last 30% of
%sessions. Columns 1 and 3: Blanco. Columns 2 and 4: Jack.
%Columns 1 and 2: Location 1 (V4). Columns 3 and 4: Location 2 (V1). 
%Rows 1, 2 & 3: proportion correct, slope, PSE.
onExternalHD=0;
plotRTerr=1;
roving=2;
plotStage3=1;
animals=[{'blanco'} {'jack'}];
if roving==0
    numconds=14;
    plotHorizontal=1;
    stage3partial=1;%0: include all 5 sessions in V4_2. 1: include just last session. 2: include last 2 sessions, etc.
    sampleContrasts=30;
elseif roving==1||roving==2 %1: normal roving task; 2: control roving task in Jack
    numconds=12;
    plotHorizontal=0;
    stage3partial=1;
    sampleContrasts=[20 30 40];
    if roving==2
        animals={'jack'};
    end
end
areas=[{'v4_1'} {'v1_1'} {'v4_2'};{'v1_2_1'} {'v1_2_2'} {'v1_2_3'};{'v1_4_1'} {'v1_4_2'} {'v1_4_3'}];
rows=4;
if plotRTerr==1
    rows=5;
end
if roving==1||roving==2
    bjboxplots_s1ls2l=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.1, 0.5, 0.8]);%[left, bottom, width, height]
    hold all
end
for sampleContrastInd=1:length(sampleContrasts)
    sampleContrast=sampleContrasts(sampleContrastInd);
    bjboxplots{sampleContrastInd}=figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.2, 0.1, 0.5, 0.8]);%[left, bottom, width, height]
    hold all
    for animalInd=1:length(animals)
        if roving==2
            monkeyMultiple=1;
        else
            monkeyMultiple=animalInd-1;%blanco:0; jack: 1
        end
        animal=animals{animalInd};
        %V4_1/V1_2_1
        area=areas{roving+1,1};
        if onExternalHD==1
            rootFolder='G:\PL_backup_060413';
        else
            rootFolder='F:';
        end
        if psychoOnly==1
            loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_',area,'_',num2str(sampleContrast),'_psycho_only.mat allMeanPerf'];
        else
            loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_',area,'_',num2str(sampleContrast),'.mat allMeanPerf'];
        end
        eval(loadText)
        maxSess=size(allMeanPerf,1);
        if strcmp(area,'v4_1')
            numN=floor((maxSess-1)*0.3);%exclude horizontal grating
        else
            numN=floor(maxSess*0.3);%number of sessions in each set
        end
        % boxplot([allMeanPerf(1:numN,2) allMeanPerf(end-numN+1:end,2)]);hold on
        grouping=[zeros(numN,1);ones(numN,1)];%for groups in boxplot
        bpc=[allMeanPerf(1:numN,2);allMeanPerf(end-numN:end-1,2)];%exclude session with horizontal Gabor
        bsl=[allMeanPerf(1:numN,4+numconds);allMeanPerf(end-numN:end-1,4+numconds)];
        bpse=[allMeanPerf(1:numN,3+numconds);allMeanPerf(end-numN:end-1,3+numconds)];
        brt=[allMeanPerf(1:numN,7+numconds)/1000;allMeanPerf(end-numN:end-1,7+numconds)/1000];
        brterr=[allMeanPerf(1:numN,23+numconds*3)/1000;allMeanPerf(end-numN:end-1,23+numconds*3)/1000];
        %stage 1 early vs late
        %     p1el(1+monkeyMultiple)=anovan(bpc(1:end),grouping(1:end));
        %     p1el(5+monkeyMultiple)=anovan(bsl(1:end),grouping(1:end));
        %     p1el(9+monkeyMultiple)=anovan(bpse(1:end),grouping(1:end));
        %     p1el(13+monkeyMultiple)=anovan(brt(1:end),grouping(1:end));
        %     p1el(17+monkeyMultiple)=anovan(brterr(1:end),grouping(1:end-1));
        pcs1h=allMeanPerf(end,2);bpcs1l=allMeanPerf(end-numN:end-1,2);%exclude session with horizontal Gabor
        sls1h=allMeanPerf(end,4+numconds);bsls1l=allMeanPerf(end-numN:end-1,4+numconds);
        pses1h=allMeanPerf(end,3+numconds);bpses1l=allMeanPerf(end-numN:end-1,3+numconds);
        rts1h=allMeanPerf(end,7+numconds)/1000;brts1l=allMeanPerf(end-numN:end-1,7+numconds)/1000;
        rterrs1h=allMeanPerf(end,23+numconds*3)/1000;brterrs1l=allMeanPerf(end-numN:end-1,23+numconds*3)/1000;
        [p,t,stats,terms]=anovan(bpc(1:end),grouping(1:end));
        p1elAll(1,1+monkeyMultiple*8:8+monkeyMultiple*8)=[mean(bpc(grouping==0)*100),var(bpc(grouping==0)*100),mean(bpcs1l*100),var(bpcs1l*100),t{2,3},t{3,3},t{2,6},p];%exclude session with horizontal Gabor
        [p,t,stats,terms]=anovan(bsl(1:end),grouping(1:end));
        p1elAll(2,1+monkeyMultiple*8:8+monkeyMultiple*8)=[mean(bsl(grouping==0)),var(bsl(grouping==0)),mean(bsls1l),var(bsls1l),t{2,3},t{3,3},t{2,6},p];
        [p,t,stats,terms]=anovan(bpse(1:end),grouping(1:end));
        p1elAll(3,1+monkeyMultiple*8:8+monkeyMultiple*8)=[mean(bpse(grouping==0)),var(bpse(grouping==0)),mean(bpses1l),var(bpses1l),t{2,3},t{3,3},t{2,6},p];
        [p,t,stats,terms]=anovan(brt(1:end),grouping(1:end));
        p1elAll(4,1+monkeyMultiple*8:8+monkeyMultiple*8)=[mean(brt(grouping==0)),var(brt(grouping==0)),mean(brts1l),var(brts1l),t{2,3},t{3,3},t{2,6},p];
        [p,t,stats,terms]=anovan(brterr(1:end),grouping(1:end));
        p1elAll(5,1+monkeyMultiple*8:8+monkeyMultiple*8)=[mean(brterr(grouping==0)),var(brterr(grouping==0)),mean(brterrs1l),var(brterrs1l),t{2,3},t{3,3},t{2,6},p];
        [h p1e1(1+monkeyMultiple) CI stats1e1(1+monkeyMultiple)]=ttest(bpc(1:numN),bpc(end-numN+1:end));
        [h p1e1(1+monkeyMultiple) CI stats1e1(3+monkeyMultiple)]=ttest(bsl(1:numN),bsl(end-numN+1:end));
        [h p1e1(1+monkeyMultiple) CI stats1e1(5+monkeyMultiple)]=ttest(bpse(1:numN),bpse(end-numN+1:end));
        [h p1e1(1+monkeyMultiple) CI stats1e1(7+monkeyMultiple)]=ttest(brt(1:numN),brt(end-numN+1:end));
        [h p1e1(1+monkeyMultiple) CI stats1e1(9+monkeyMultiple)]=ttest(brterr(1:numN),brterr(end-numN+1:end));
        r1lh{1+monkeyMultiple}=[min(bpcs1l),max(bpcs1l),pcs1h];
        r1lh{3+monkeyMultiple}=[min(bsls1l),max(bsls1l),sls1h];
        r1lh{5+monkeyMultiple}=[min(bpses1l),max(bpses1l),pses1h];
        r1lh{7+monkeyMultiple}=[min(brts1l),max(brts1l),rts1h];
        r1lh{9+monkeyMultiple}=[min(brterrs1l),max(brterrs1l),rterrs1h];
        %range of stage 1 late vs horizontal stimulus session
        [h p1lh(1+monkeyMultiple)]=ttest(bpcs1l,pcs1h);
        [h p1lh(3+monkeyMultiple)]=ttest(bsls1l,sls1h);
        [h p1lh(5+monkeyMultiple)]=ttest(bpses1l,pses1h);
        [h p1lh(7+monkeyMultiple)]=ttest(brts1l,rts1h);
        [h p1lh(9+monkeyMultiple)]=ttest(brterrs1l,rterrs1h);
        
        if plotStage3==1%V4_2/V1_2_3
            area=areas{roving+1,3};
            loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_',area,'_',num2str(sampleContrast),'.mat allMeanPerf'];
            eval(loadText)
            maxSess=size(allMeanPerf,1);
            if stage3partial==0
                grouping=[grouping;zeros(maxSess,1)+2];
                bpc=[bpc;allMeanPerf(:,2)];
                bsl=[bsl;allMeanPerf(:,4+numconds)];
                bpse=[bpse;allMeanPerf(:,3+numconds)];
                brt=[brt;allMeanPerf(:,7+numconds)/1000];
                brterr=[brterr;allMeanPerf(:,23+numconds*3)/1000];
            else
                grouping=[grouping;zeros(stage3partial,1)+2];
                bpc=[bpc;allMeanPerf(end-stage3partial+1:end,2)];
                bsl=[bsl;allMeanPerf(end-stage3partial+1:end,4+numconds)];
                bpse=[bpse;allMeanPerf(end-stage3partial+1:end,3+numconds)];
                brt=[brt;allMeanPerf(end-stage3partial+1:end,7+numconds)/1000];
                brterr=[brterr;allMeanPerf(end-stage3partial+1:end,23+numconds*3)/1000];
                if stage3partial==1
                    finalSessionVals{sampleContrastInd}(:,monkeyMultiple+1)=[bpc(end);bsl(end);bpse(end);brt(end);brterr(end)];
                end
            end
            %stage 1 late vs stage 3
            p1l3(1+monkeyMultiple)=anovan(bpc(grouping>0),grouping(grouping>0));%compare data with grouping values of 1 and 2
            p1l3(3+monkeyMultiple)=anovan(bsl(grouping>0),grouping(grouping>0));
            p1l3(5+monkeyMultiple)=anovan(bpse(grouping>0),grouping(grouping>0));
            p1l3(7+monkeyMultiple)=anovan(brt(grouping>0),grouping(grouping>0));
            p1l3(9+monkeyMultiple)=anovan(brterr(grouping>0),grouping(grouping>0));
            meanperf3(1,monkeyMultiple+1)=mean(allMeanPerf(end-stage3partial+1:end,2));
            meanperf3(2,monkeyMultiple+1)=mean(allMeanPerf(end-stage3partial+1:end,4+numconds));
            meanperf3(3,monkeyMultiple+1)=mean(allMeanPerf(end-stage3partial+1:end,3+numconds));
            meanperf3(4,monkeyMultiple+1)=mean(allMeanPerf(end-stage3partial+1:end,7+numconds)/1000);
            meanperf3(5,monkeyMultiple+1)=mean(allMeanPerf(end-stage3partial+1:end,23+numconds*3)/1000);
        end
        figure(bjboxplots{sampleContrastInd})
        subplot(rows,4,1+monkeyMultiple);%3: pc, slope, PSE. 4: 2 monkeys and 2 locations
        boxplot(bpc(1:end-1),grouping(1:end-1));hold on
        %     text('String','x','Position',[2.25,bpc(end)],'Color','k');
        if plotHorizontal==1
            plot(2.3,pcs1h,'Marker','o');
        end
        %     xlabel('sessions');
        %     ylabel('Proportion correct');
        %     set(gca,'XTick',[1;2])
        %     set(gca,'XTickLabel',{'early';'late'})
        %     title('Stage 1 Proportion correct');
        [yAxis]=get(gca,'yLim');
        set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
        %     set(gca,'XTick',[1;2;3])
        %     set(gca,'XTickLabel',{'1 early';'1 late';'3'})
        %     title('Stages 1 and 3 proportion correct');
        
        subplot(rows,4,5+monkeyMultiple);
        boxplot(bsl(1:end-1),grouping(1:end-1));hold on
        %     text('String','x','Position',[2.25,bsl(end)],'Color','k');
        if plotHorizontal==1
            plot(2.3,sls1h,'Marker','o');
        end
        %     xlabel('Sessions');
        %     ylabel('Slope');
        %     set(gca,'XTick',[1;2])
        %     set(gca,'XTickLabel',{'early';'late'})
        %     title('Stage 1 Slope');
        [yAxis]=get(gca,'yLim');
        set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
        %     set(gca,'XTick',[1;2;3])
        %     set(gca,'XTickLabel',{'1 early';'1 late';'3'})
        %     title('Stages 1 and 3 slope');
        
        subplot(rows,4,9+monkeyMultiple);
        boxplot(bpse(1:end-1),grouping(1:end-1));hold on
        %     text('String','x','Position',[2.25,bpse(end)],'Color','k');
        if plotHorizontal==1
            plot(2.3,pses1h,'Marker','o');
        end
        %     xlabel('Sessions');
        %     ylabel('PSE');
        %     set(gca,'XTick',[1;2])
        %     set(gca,'XTickLabel',{'early';'late'})
        %     title('Stage 1 PSE');
        [yAxis]=get(gca,'yLim');
        set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
        %     set(gca,'XTick',[1;2;3])
        %     set(gca,'XTickLabel',{'1 early';'1 late';'3'})
        %     title('Stages 1 and 3 PSE');
        
        subplot(rows,4,13+monkeyMultiple);
        boxplot(brt(1:end-1),grouping(1:end-1));hold on
        %     text('String','x','Position',[2.25,brt(end)],'Color','k');
        if plotHorizontal==1
            plot(2.3,rts1h,'Marker','o');
        end
        %     xlabel('Sessions');
        %     ylabel('RT correct');
        %     set(gca,'XTick',[1;2])
        %     set(gca,'XTickLabel',{'early';'late'})
        %     title('Stage 1 RT correct');
        [yAxis]=get(gca,'yLim');
        set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
        %     set(gca,'XTick',[1;2;3])
        %     set(gca,'XTickLabel',{'1 early';'1 late';'3'})
        %     title('Stages 1 and 3 RT');
        
        if plotRTerr==1
            subplot(rows,4,17+monkeyMultiple);
            boxplot(brterr(1:end-1),grouping(1:end-1));hold on
            %     text('String','x','Position',[2.25,brt(end)],'Color','k');
            if plotHorizontal==1
                plot(2.3,rterrs1h,'Marker','o');
            end
            %         xlabel('Sessions');
            %         ylabel('RT error');
            %         set(gca,'XTick',[1;2])
            %         set(gca,'XTickLabel',{'early';'late'})
            %         title('Stage 1 RT error');
            [yAxis]=get(gca,'yLim');
            set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
            %     set(gca,'XTick',[1;2;3])
            %     set(gca,'XTickLabel',{'1 early';'1 late';'3'})
            %     title('Stages 1 and 3 RT');
        end
        
        %V1/V1_2_2
        area=areas{roving+1,2};
        loadText=['load ',rootFolder,'\PL\psycho_data\',animal,'\allMeanPerf\allMeanPerf_',area,'_',num2str(sampleContrast),'.mat allMeanPerf'];
        eval(loadText)
        maxSess=size(allMeanPerf,1);
        numN=floor(maxSess*0.3);%number of sessions in each set
        % boxplot([allMeanPerf(1:numN,2) allMeanPerf(end-numN+1:end,2)]);hold
        % on
        grouping=[zeros(numN,1);ones(numN,1)];
        bpc=[allMeanPerf(1:numN,2);allMeanPerf(end-numN+1:end,2)];
        bsl=[allMeanPerf(1:numN,4+numconds);allMeanPerf(end-numN+1:end,4+numconds)];
        bpse=[allMeanPerf(1:numN,3+numconds);allMeanPerf(end-numN+1:end,3+numconds)];bpse=real(bpse);
        brt=[allMeanPerf(1:numN,7+numconds)/1000;allMeanPerf(end-numN+1:end,7+numconds)/1000];
        brterr=[allMeanPerf(1:numN,23+numconds*3)/1000;allMeanPerf(end-numN+1:end,23+numconds*3)/1000];
        %stage 2 early vs late
        %     p1el(3+monkeyMultiple)=anovan(bpc,grouping);
        %     p1el(7+monkeyMultiple)=anovan(bsl,grouping);
        %     p1el(11+monkeyMultiple)=anovan(bpse,grouping);
        %     p1el(15+monkeyMultiple)=anovan(brt,grouping);
        %     p1el(19+monkeyMultiple)=anovan(brterr,grouping);
        [p,t,stats,terms]=anovan(bpc,grouping);
        p1elAll(1,17+monkeyMultiple*8:24+monkeyMultiple*8)=[mean(bpc(grouping==0)*100),var(bpc(grouping==0)*100),mean(bpc(grouping==1)*100),var(bpc(grouping==1)*100),t{2,3},t{3,3},t{2,6},p];%exclude session with horizontal Gabor
        [p,t,stats,terms]=anovan(bsl,grouping);
        p1elAll(2,17+monkeyMultiple*8:24+monkeyMultiple*8)=[mean(bsl(grouping==0)),var(bsl(grouping==0)),mean(bsl(grouping==1)),var(bsl(grouping==1)),t{2,3},t{3,3},t{2,6},p];
        [p,t,stats,terms]=anovan(bpse,grouping);
        p1elAll(3,17+monkeyMultiple*8:24+monkeyMultiple*8)=[mean(bpse(grouping==0)),var(bpse(grouping==0)),mean(bpse(grouping==1)),var(bpse(grouping==1)),t{2,3},t{3,3},t{2,6},p];
        [p,t,stats,terms]=anovan(brt,grouping);
        p1elAll(4,17+monkeyMultiple*8:24+monkeyMultiple*8)=[mean(brt(grouping==0)),var(brt(grouping==0)),mean(brt(grouping==1)),var(brt(grouping==1)),t{2,3},t{3,3},t{2,6},p];
        [p,t,stats,terms]=anovan(brterr,grouping);
        p1elAll(5,17+monkeyMultiple*8:24+monkeyMultiple*8)=[mean(brterr(grouping==0)),var(brterr(grouping==0)),mean(brterr(grouping==1)),var(brterr(grouping==1)),t{2,3},t{3,3},t{2,6},p];
        [h p1e1(1+monkeyMultiple) CI stats1e1(2,1+monkeyMultiple)]=ttest(bpc(1:numN),bpc(end-numN+1:end));
        [h p1e1(1+monkeyMultiple) CI stats1e1(2,3+monkeyMultiple)]=ttest(bsl(1:numN),bsl(end-numN+1:end));
        [h p1e1(1+monkeyMultiple) CI stats1e1(2,5+monkeyMultiple)]=ttest(bpse(1:numN),bpse(end-numN+1:end));
        [h p1e1(1+monkeyMultiple) CI stats1e1(2,7+monkeyMultiple)]=ttest(brt(1:numN),brt(end-numN+1:end));
        [h p1e1(1+monkeyMultiple) CI stats1e1(2,9+monkeyMultiple)]=ttest(brterr(1:numN),brterr(end-numN+1:end));
        if roving==1||roving==2%compare flanker late with flankerless late (late Stage 4 vs late Stage 5)
            figure(bjboxplots_s1ls2l)
            groupings1ls2l=[zeros(length(bpcs1l),1);ones(length(bpc(grouping==1)),1)];
            bpcs1ls21=[bpcs1l*100;bpc(grouping==1)*100];
            bsls1ls21=[bsls1l;bsl(grouping==1)];
            bpses1ls21=[bpses1l;bpse(grouping==1)];
            brts1ls21=[brts1l;brt(grouping==1)];
            brterrs1ls21=[brterrs1l;brterr(grouping==1)];
            [p,t,stats,terms]=anovan(bpcs1ls21,groupings1ls2l);
            p1l2lAll{sampleContrastInd,monkeyMultiple+1}(1,:)=[mean(bpcs1l*100),var(bpcs1l*100),mean(bpc(grouping==1)*100),var(bpc(grouping==1)*100),t{2,3},t{3,3},t{2,6},p];%exclude session with horizontal Gabor
            [p,t,stats,terms]=anovan(bsls1ls21,groupings1ls2l);
            p1l2lAll{sampleContrastInd,monkeyMultiple+1}(2,:)=[mean(bsls1l),var(bsls1l),mean(bsl(grouping==1)),var(bsl(grouping==1)),t{2,3},t{3,3},t{2,6},p];
            [p,t,stats,terms]=anovan(bpses1ls21,groupings1ls2l);
            p1l2lAll{sampleContrastInd,monkeyMultiple+1}(3,:)=[mean(bpses1l),var(bpses1l),mean(bpse(grouping==1)),var(bpse(grouping==1)),t{2,3},t{3,3},t{2,6},p];
            [p,t,stats,terms]=anovan(brts1ls21,groupings1ls2l);
            p1l2lAll{sampleContrastInd,monkeyMultiple+1}(4,:)=[mean(brts1l),var(brts1l),mean(brt(grouping==1)),var(brt(grouping==1)),t{2,3},t{3,3},t{2,6},p];
            [p,t,stats,terms]=anovan(brterrs1ls21,groupings1ls2l);
            p1l2lAll{sampleContrastInd,monkeyMultiple+1}(5,:)=[mean(brterrs1l),var(brterrs1l),mean(brterr(grouping==1)),var(brterr(grouping==1)),t{2,3},t{3,3},t{2,6},p];
            subplot(rows,length(sampleContrasts),sampleContrastInd);
            boxplot(bpcs1ls21,groupings1ls2l);
            [yAxis]=get(gca,'yLim');
            set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
            
            subplot(rows,length(sampleContrasts),sampleContrastInd+3);
            boxplot(bsls1ls21,groupings1ls2l);
            [yAxis]=get(gca,'yLim');
            set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
            
            subplot(rows,length(sampleContrasts),sampleContrastInd+6);
            boxplot(bpses1ls21,groupings1ls2l);
            [yAxis]=get(gca,'yLim');
            set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
            
            subplot(rows,length(sampleContrasts),sampleContrastInd+9);
            boxplot(brts1ls21,groupings1ls2l);
            [yAxis]=get(gca,'yLim');
            set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
            
            if plotRTerr==1
                subplot(rows,length(sampleContrasts),sampleContrastInd+12);
                boxplot(brterrs1ls21,groupings1ls2l);
                [yAxis]=get(gca,'yLim');
                set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
            end
        end
        figure(bjboxplots{sampleContrastInd})
        subplot(rows,4,3+monkeyMultiple);
        boxplot(bpc,grouping);
        [yAxis]=get(gca,'yLim');
        set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
        
        subplot(rows,4,7+monkeyMultiple);
        boxplot(bsl,grouping);
        [yAxis]=get(gca,'yLim');
        set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
        
        subplot(rows,4,11+monkeyMultiple);
        boxplot(bpse,grouping);
        [yAxis]=get(gca,'yLim');
        set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
        
        subplot(rows,4,15+monkeyMultiple);
        boxplot(brt,grouping);
        [yAxis]=get(gca,'yLim');
        set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
        
        if plotRTerr==1
            subplot(rows,4,19+monkeyMultiple);
            boxplot(brterr,grouping);
            %         xlabel('Sessions');
            %         ylabel('RT error');
            %         set(gca,'XTick',[1;2])
            %         set(gca,'XTickLabel',{'2 early';'2 late'})
            %         title('Stage 2 RT error');
            [yAxis]=get(gca,'yLim');
            set(gca,'yLim',[yAxis(1) (yAxis(2)-yAxis(1))*1.2+yAxis(1)]);
        end
    end
    p1elAllSamples{sampleContrastInd}=p1elAll;
    %p1el=reshape(p1el,4,5)'
    p1l3samples{sampleContrastInd}=reshape(p1l3,2,5)'
    if plotStage3==1
        r1lhtemp=reshape(r1lh,2,5)'
        p1lh=reshape(p1lh,2,5)'
        if roving==0
            r1lhFormatted=[];
            for i=1:5
                r1lhFormatted=[r1lhFormatted;r1lhtemp{i,1} r1lhtemp{i,2}];
            end
        elseif roving==1
            r1lhFormattedTemp=[];
            for i=1:5
                r1lhFormattedTemp=[r1lhFormattedTemp;r1lhtemp{i,1} r1lhtemp{i,2}];
            end
            r1lhFormatted{sampleContrastInd}=r1lhFormattedTemp;
        elseif roving==2
            r1lhFormattedTemp=[];
            for i=1:5
                r1lhFormattedTemp=[r1lhFormattedTemp;r1lhtemp{i,2}];
            r1lhFormatted{sampleContrastInd}=r1lhFormattedTemp;
            end
        end
    end
end
close all hidden
