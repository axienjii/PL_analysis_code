function read_bj_SE_xcorr6skew2(animal,area)
%Written by Xing 16/09/10, slightly modified on 20/09/10 to accept 2nd
%input arg.
%Modified from read_blanco_SE_xcorr5skew.
%Reads values of proportions of corr coefs which lie within 95% or 99%
%CIs of the bootstrapped data calculated for each session, for both the
%within-cell-between-session coefs, and the across-cell-and-session coefs.
%Plots proportions for each session in a subplot within a big diagram which
%displays all sessions together. Saves 2 images to file, e.g.
%95_CI_barplots.
%Also generates graph for each session, indicating whether the coefs
%which compare that session with the others fall within the CIs for the
%bootstrapped data from that session. Draws a red cross if within limits,
%and a grey cross if not. Saves image to file, e.g. 8_CI_table.

minusSpontan=0;
sigma=8;
if minusSpontan==1
    subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
elseif minusSpontan==0
    subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_wspontan_',area];%folder for stimulus-evoked responses without any subtraction of spontaneous activity levels
end

matTexts=[{'trans_allDist'} {'grandTrialsList'} {'CItable'}];
for i=1:length(matTexts)
    matName=[matTexts{i},'_',area];
    matPath=fullfile('F:','PL','xcorr',animal,subfolder,matName);
    loadText=['load ',matPath,' ',matTexts{i}];
    eval(loadText);
end
allDist=trans_allDist;

sessionNums=main_raw_sessions_final(animal,area,[],0);
channels=[grandTrialsList{:,1}];
ind=[];
% for i=1:size(CItable,1)
for i=1:size(CItable,1)
    if ~isempty(find(CItable(i,:,1)~=-1,1))
        ind=[ind i];
    end
end
newCItable=CItable((ind),:,:);
channels=channels(ind);
%get rid of sessions for which no data is present:
ind=[];
for i=1:size(CItable,2)
    if ~isempty(find(CItable(:,i,1)~=-1,1))
        ind=[ind i];
    end
end
newCItable=newCItable(:,(ind),:);
sessionNums=sessionNums(ind);
printFlag=0;
figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05,0.05, 0.9, 0.9]);
for i=1:size(newCItable,1)
    channel=channels(i);
    if i>floor(size(newCItable,1)/2)
        rowIndex=i-floor(size(newCItable,1)/2);
    else
        rowIndex=i;
    end
    if i==floor(size(newCItable,1)/2)+1
        if printFlag==0
            printFlag=1;
            set(gcf,'Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05,0.05, 0.9, 0.9]);
            ptext=sprintf('w: comparing sessions within a cell         b: across cells');
            text('Position',[-15 -2],'FontSize',9,'String',ptext);
            imageName=['95_CI_barplots1_',area];
            imageFolder=fullfile('F:','PL','xcorr',animal,subfolder);
            imagePath=fullfile(imageFolder,imageName);
            printtext=['print -dpng ',imagePath];
            set(gcf,'PaperPositionMode','auto')
            eval(printtext);
            figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05,0.05, 0.9, 0.9]);
        end
    end
    for j=1:size(newCItable,2)        
        subplot(size(newCItable,1)-floor(size(newCItable,1)/2),size(newCItable,2),j+size(newCItable,2)*(rowIndex-1));
        if newCItable(i,j,1)==-1
            axis off
        end
        if j==1
            if newCItable(i,j,1)==-1
                text('Position',[-1 0.5],'FontSize',9,'String',num2str(channel));hold on
            end
            if newCItable(i,j,1)~=-1
                text('Position',[-4 0.5],'FontSize',9,'String',num2str(channel));hold on
            end
        end
        if newCItable(i,j,1)~=-1
            proportions=[newCItable(i,j,1) newCItable(i,j,2)];
            bar([0 1],proportions);
            set(gca,'YLim',[0 1]);
            set(gca,'YTick',[0 1]);
            set(gca,'YTickLabel',[0 1]);
            set(gca,'XTick',[0 1]);
            set(gca,'XLim',[-1 2]);
            set(gca,'XTickLabel',['w';'b']);
            if i==1||i==floor(size(newCItable,1)/2)+1
                title(sessionNums(j));
                if j==1
                    ptext=sprintf('proportion of correlation coefs within CI     %s',matPath);
                    text('Position',[-2 2.5],'FontSize',9,'String',ptext);
                end
            end
        end
    end
end
ptext=sprintf('w: comparing sessions within a cell         b: across cells');
text('Position',[-15 -2],'FontSize',9,'String',ptext);
imageName=['95_CI_barplots2_',area];
imageFolder=fullfile('F:','PL','xcorr',animal,subfolder);
imagePath=fullfile(imageFolder,imageName);
printtext=['print -dpng ',imagePath];
set(gcf,'PaperPositionMode','auto')
eval(printtext);
close all

printFlag=0;
figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05,0.05, 0.9, 0.9]);
for i=1:size(newCItable,1)
    channel=channels(i);
    if i>floor(size(newCItable,1)/2)
        rowIndex=i-floor(size(newCItable,1)/2);
    else
        rowIndex=i;
    end
    if i==floor(size(newCItable,1)/2)+1
        if printFlag==0
            printFlag=1;
            set(gcf,'Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05,0.05, 0.9, 0.9]);
            ptext=sprintf('w: comparing sessions within a cell         b: across cells');
            text('Position',[-15 -5],'FontSize',9,'String',ptext);
            imageName=['99_CI_barplots1_',area];
            imageFolder=fullfile('F:','PL','xcorr',animal,subfolder);
            imagePath=fullfile(imageFolder,imageName);
            printtext=['print -dpng ',imagePath];
            set(gcf,'PaperPositionMode','auto')
            eval(printtext);
            figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05,0.05, 0.9, 0.9]);
        end
    end
    for j=1:size(newCItable,2)
        subplot(size(newCItable,1)-floor(size(newCItable,1)/2),size(newCItable,2),j+size(newCItable,2)*(rowIndex-1));
        if newCItable(i,j,1)==-1
            axis off
        end
        if j==1
            if newCItable(i,j,1)==-1
                text('Position',[-1 0.5],'FontSize',9,'String',num2str(channel));hold on
            end
            if newCItable(i,j,1)~=-1
                text('Position',[-4 0.5],'FontSize',9,'String',num2str(channel));hold on
            end
        end
        if newCItable(i,j,1)~=-1
            proportions=[newCItable(i,j,3) newCItable(i,j,4)];
            bar([0 1],proportions);
            set(gca,'YLim',[0 1]);
            set(gca,'YTick',[0 1]);
            set(gca,'YTickLabel',[0 1]);
            set(gca,'XTick',[0 1]);
            set(gca,'XLim',[-1 2]);
            set(gca,'XTickLabel',['w';'b']);
            if i==1||i==floor(size(newCItable,1)/2)+1
                title(sessionNums(j));
                if j==1
                    ptext=sprintf('proportion of correlation coefs within CI     %s',matPath);
                    text('Position',[-2 3.5],'FontSize',9,'String',ptext);
                end
            end
        end
    end
end
ptext=sprintf('w: comparing sessions within a cell         b: across cells');
text('Position',[-15 -5],'FontSize',9,'String',ptext);
imageName=['99_CI_barplots2_',area];
imageFolder=fullfile('F:','PL','xcorr',animal,subfolder);
imagePath=fullfile(imageFolder,imageName);
printtext=['print -dpng ',imagePath];
set(gcf,'PaperPositionMode','auto')
eval(printtext);
close all

numStds=[1.96 2.576];percent=[95 99];
for stdCount=1:length(numStds)
    for i=1:size([grandTrialsList{:,1}],2)%each cell
        figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.1,0.3, 0.5, 0.4]);
        if ~isempty(grandTrialsList{i,2})
            yCount=1;ySessions=[];
            for rowCount=1:size(allDist,1)
                if ~isempty(allDist{rowCount,4})
                    if allDist{rowCount,1}==grandTrialsList{i,1}
                        meanBootC=mean(allDist{rowCount,3});
                        stdC=std(allDist{rowCount,3});
                        all_sessions=grandTrialsList{i,3};
                        all_sessions(all_sessions==3551)=355;all_sessions(all_sessions==3552)=355.5;
                        orderedSessions=1:length(all_sessions);
                        sessions=orderedSessions(all_sessions~=allDist{rowCount,2});
                        withinCount=0;
                        for k=1:length(allDist{rowCount,4})
                            if allDist{rowCount,4}(k)<=meanBootC+numStds(stdCount)*stdC&&allDist{rowCount,4}(k)>=meanBootC-numStds(stdCount)*stdC%check whether value of coef falls within 2 SDs
                                plot(sessions(k),yCount,'rx');hold on
                                withinCount=withinCount+1;
                            else
                                plot(sessions(k),yCount,'x','Color',[0.7 0.7 0.7]);hold on
                            end
                        end
                        text('Position',[length(all_sessions)+1+0.03*(length(all_sessions)+2) yCount],'FontSize',8,'String',num2str(withinCount),'Color','k');
                        yCount=yCount+1;
                        ySessions=[ySessions orderedSessions(find(all_sessions==allDist{rowCount,2}))];
                    end
                end
            end
            ylabel('session for which CI is calculated');xlabel('compared with session');
            set(gca,'YTick',1:1:yCount-1);
            set(gca,'YTickLabel',ySessions);
            set(gca,'YLim',[0 yCount]);
            set(gca,'YDir','reverse')
            set(gca,'XTick',1:length(all_sessions));
            set(gca,'XTickLabel',1:length(all_sessions));
            set(gca,'XLim',[0 length(all_sessions)+1]);
            ptext=sprintf('%s            # of trials per bootstrapped CI: %d',num2str(grandTrialsList{i,1}),length(allDist{i,3}));
            text('Position',[-1 -(yCount-1)/20],'FontSize',9,'String',ptext);
            ptext2=sprintf('outside %d%% CI',percent(stdCount));
            ptext3=sprintf('within   %d%% CI',percent(stdCount));
            text('Position',[length(all_sessions)/2 -(yCount-1)/40],'FontSize',9,'String',ptext2,'Color','k');
            text('Position',[length(all_sessions)/2 -(yCount-1)/15],'FontSize',9,'String',ptext3,'Color','r');
            tallyText=sprintf('tally (%d)',length(all_sessions)-1);
            text('Position',[length(all_sessions)+1+0.02*(length(all_sessions)+2) 0],'FontSize',8,'String',tallyText,'Color','k');
            channel=num2str(grandTrialsList{i,1});
            imageName=[num2str(channel),'_',num2str(percent(stdCount)),'_',area,'_CI_table'];
            imageFolder=fullfile('F:','PL','xcorr',animal,subfolder);
            imagePath=fullfile(imageFolder,imageName);
            printtext=['print -dpng ',imagePath];
            set(gcf,'PaperPositionMode','auto')
            eval(printtext);
            close all
        end
    end
end
