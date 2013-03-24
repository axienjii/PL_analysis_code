function read_bj_SE_xcorr7skew(animal,area)
%Written by Xing 16/09/10, slightly modified on 20/09/10 to accept 2nd
%input arg.
%Modified from read_blanco_SE_xcorr6skew2.
%Reads values of proportions of corr coefs which lie within 95% or 99%
%CIs of the bootstrapped data calculated for each session, for both the
%within-cell-between-session coefs, and the across-cell-and-session coefs.
%Generates plots proportions with within-cell values on the y-axis, and between-cell values on the x-axis.
%Saves image to file.

minusSpontan=0;
sigma=8;

intervalSizes=[95 99];
for intervalInd=1:length(intervalSizes)
    figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05,0.05, 0.9, 0.9]);
    set(gcf,'Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05,0.05, 0.3, 0.3]);
    animalTexts=[{'subject B'} {'subject J'}];
    animals=[{'blanco'} {'jack'}];
    areaTexts=[{'V4'} {'V1'}];
    areas=[{'v4_1'} {'v1_1'}];
    for animalInd=1:2
        animal=animals{animalInd};
        for areaInd=1:2
            area=areas{areaInd};
            folder=fullfile('F:','PL','xcorr',animal);
            if minusSpontan==1
                subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
            elseif minusSpontan==0
                subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_wspontan_',area];%folder for stimulus-evoked responses without any subtraction of spontaneous activity levels
            end
            
            matTexts=[{'trans_allDist'} {'grandTrialsList'} {'CItable'}];
            for i=1:length(matTexts)
                matName=[matTexts{i},'_',area];
                matPath=fullfile(folder,subfolder,matName);
                loadText=['load ',matPath,' ',matTexts{i}];
                eval(loadText);
            end
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
            count=0;
            higherCount=0;
            subplotInd=subplot(2,2,animalInd+2*(areaInd-1));
            for i=1:size(newCItable,1)
                channel=channels(i);
                ptext=sprintf('w: comparing sessions within a cell         b: across cells');
                text('Position',[-15 -2],'FontSize',9,'String',ptext);
                for j=1:size(newCItable,2)
                    %         subplot(size(newCItable,1)-floor(size(newCItable,1)/2),size(newCItable,2),j+size(newCItable,2)*(rowIndex-1));
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
                        proportions=[newCItable(i,j,1+(intervalInd-1)*2) newCItable(i,j,2+(intervalInd-1)*2)];
                        plot(proportions(2),proportions(1),'ko','MarkerSize',3);hold on
                        count=count+1;
                        if proportions(2)<proportions(1)
                            higherCount=higherCount+1;
                        end
                    end
                end
            end
            set(gca,'YLim',[0 1]);
            set(gca,'YTick',[0 1]);
            set(gca,'YTickLabel',[0 1]);
            set(gca,'XTick',[0 1]);
            set(gca,'XLim',[0 1]);
            set(gca,'XTickLabel',[0 1]);
            line([0 1],[0 1],'LineStyle',':','Color','k','LineWidth',1);
            axis square
            if animalInd==1
                ptext=areaTexts{areaInd};
                text('Position',[-0.5 0.5],'FontSize',9,'String',ptext);
                ylabel('proportion of R_a within CI');
            end
            if areaInd==1
                title(animalTexts{animalInd});
            else
                xlabel('proportion of R_c within CI');
            end
            tallyHigherWithin(intervalInd,animalInd+2*(areaInd-1))={[higherCount count]};%number of cases where proportion within-cell R-values are higher than between-cell R-values
        end
    end
    % ptext=sprintf('w: comparing sessions within a cell         b: across cells');
    % text('Position',[-15 -2],'FontSize',9,'String',ptext);
    imageName=[num2str(intervalSizes(intervalInd)),'_CI_scatterplots_V4_V1'];
    imageFolder=fullfile('F:','PL','xcorr');
    imagePath=fullfile(imageFolder,imageName);
    printtext=['print -dpng -r300 ',imagePath];
    set(gcf,'PaperPositionMode','auto')
    eval(printtext);
end

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
            imageName=['99_CI_scatterplots1_',area];
            imageFolder=fullfile('F:','PL','xcorr',animal,subfolder);
            imagePath=fullfile(imageFolder,imageName);
            printtext=['print -dpng ',imagePath];
            set(gcf,'PaperPositionMode','auto')
            eval(printtext);
            figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.05,0.05, 0.9, 0.9]);
        end
    end
    for j=1:size(newCItable,2)
        %         subplot(size(newCItable,1)-floor(size(newCItable,1)/2),size(newCItable,2),j+size(newCItable,2)*(rowIndex-1));
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
            bar(proportions(2),proportions(2));
            set(gca,'YLim',[0 1]);
            set(gca,'YTick',[0 1]);
            set(gca,'YTickLabel',[0 1]);
            set(gca,'XTick',[0 1]);
            set(gca,'XLim',[0 1]);
            set(gca,'XTickLabel',[0 1]);
            xlabel('R_c');
            ylabel('R_a');
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