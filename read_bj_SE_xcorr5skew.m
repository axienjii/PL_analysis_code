function read_bj_SE_xcorr5skew(animal,area,sessionNums)
%Written by Xing 16/09/10, slightly modified on 20/09/10 to accept 2nd
%input arg. Further modified on 18/10/10 to calculate best fit line to data
%points with polyfit and write slope to 7th column in allDist.mat array,
%and intercept to 8th column. Done for each session for each cell. Slopes examined later by
%read_blanco_SE_xcorr5skewslope for signs of significant changes over time,
%which would imply changes in amount of correlation occur over time. 
%Modified from read_blanco_SE_xcorr5.
%Transforms data (corr coef values) because the original coef values have a
%negatively (left) skewed distribution, and not possible to calculate CIs
%using non-normally-distributed data.
%Tried various types of transformations and tested for normality with
%lillietest function, also calculated skewness and kurtosis. Found that the
%square root transform does pretty well (after data have undergone a
%preliminary transformation to make them positively skewed), and reduces
%skewness substantially. Performs rest of analysis as does
%read_blanco_SE_xcorr5.
%Reads values of correlation coefficients for several sets of data:
%bootstrapped from within-session for each cell, real correlation coefs
%between sessions for each cells, real coefs across cells and sessions, and
%if desired- more for curiousity's sake than anything else- real coefs
%across cells within session.
%Plots CIs (95%) and mean of bootstrapped coefs on graph and draws data
%points for correlation coef value for each comparison between PSTH values for that session
%and the other sessions for that cell. Also draws mean and CIs of
%across-cell comparisons, which might not be quite the right thing to do-
%would it be better to generate bootstrapped values for this data set? Not
%sure! Saves image to file, e.k. 8_343_CI_PSTHs, in folder PSTH_xcorr_sm6ms.
%Modified from read_blanco_SE_xcorr4 which performed analysis for each good
%session of each cell.
%Reads .mat array allDist.mat.
%Ch # in 1st column, session # in 2nd, and the 4 arrays of coef values in 'cell' format
%in the 3rd to 6th columns: 3. coefs from
%bootstrapped  values for that session, coefs from comparisons between: 4.
%within-cell-between-sessions, 5. across cells across sessions, 6. across
%cells within session. Data compiled across good sessions from top to
%bottom row of saved variable 'allDist.'

minusSpontan=0;
sigma=8;
if minusSpontan==1
    subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
elseif minusSpontan==0
    subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_wspontan_',area];%folder for stimulus-evoked responses without any subtraction of spontaneous activity levels
end
matGrandName=['grandTrialsList_',area];
matGrandPath=fullfile('F:','PL','xcorr',animal,subfolder,matGrandName);
loadText=['load ',matGrandPath,' grandTrialsList'];
eval(loadText);

matDistName=['trans_allDist_',area];
matDistPath=fullfile('F:','PL','xcorr',animal,subfolder,matDistName);
loadText=['load ',matDistPath,' trans_allDist'];
eval(loadText);
allDist=trans_allDist;

if nargin<3 || isempty(sessionNums)
    sessionNums = main_raw_sessions_final(animal,area,[],0);
end
numSessions=length(sessionNums);
CItable=zeros(size(grandTrialsList,1),length(sessionNums),4)-1;
numStds=[1.96];
stdDevText=sprintf('- - std dev %s %s %s',num2str(numStds));
skewedBoot=[];
for i=1:size(allDist,1)
    if ~isempty(allDist{i,4})
%transform data from skewed distribution to one with normal distribution
%(distribution of coefs is negatively skewed (left skewed)
% transformedBootDist=1-allDist{i,3};%another attempt- do a reversal by subtracting each value from highest value, setting range from 0 to 1
% transformedBootDist=sqrt(transformedBootDist);%then square the data
% transformedBootDist=1-transformedBootDist;%perform reversal to restore correct order of values
% % % transformedBootDist=max(allDist{i,3})-allDist{i,3};%another attempt- do a reversal by subtracting each value from highest value, sets range relative to max coef value (somewhat differently from above)
% % % transformedBootDist=sqrt(transformedBootDist);%then square the data
% % %less successful transformation methods:
% % %         transformedBootDist=allDist{i,3}/median(allDist{i,3});%attempt at
% % %         transforming data with skewed distribution to normal distribution
% % %         transformedBootDist=(allDist{i,3}).^2;%another attempt
% % %do a reversal by subtracting each value from highest value, then do square root transform:
% % %transformedBootDist=max(allDist{i,3})-allDist{i,3};
% % % transformedBootDist=nonzeros(transformedBootDist);
% % % transformedBootDist=1./transformedBootDist;
        [skew p]=lillietest(allDist{i,3});
        howSkewed=skewness(allDist{i,3});
        k=kurtosis(allDist{i,3});
%         [skew p]=lillietest(transformedBootDist);
%         howSkewed=skewness(transformedBootDist);
%         k=kurtosis(transformedBootDist);
        skewedBoot=[skewedBoot;allDist{i,1} allDist{i,2} skew p howSkewed k];
        
        figure('Color',[1,1,1],'Units', 'Normalized', 'Position',[0.1,0.3, 0.5, 0.4]);
%         meanBootC=mean(transformedBootDist);
%         stdC=std(transformedBootDist);
        meanBootC=mean(allDist{i,3});
        stdC=std(allDist{i,3});
        rowInd=find([grandTrialsList{:,1}]==allDist{i,1});
        all_sessions=grandTrialsList{rowInd,3};
        sessions=all_sessions(all_sessions~=allDist{i,2});
        all_sessions(all_sessions==3551)=355;all_sessions(all_sessions==3552)=355.5;
        sessions(sessions==3551)=355;sessions(sessions==3552)=355.5;
        line([0 numSessions],[meanBootC meanBootC],'LineStyle','-','Color','k');hold on
        for k=1:length(numStds)
            line([0 numSessions],[meanBootC+numStds(k)*stdC meanBootC+numStds(k)*stdC],'LineStyle','--','Color','k');hold on
            line([0 numSessions],[meanBootC-numStds(k)*stdC meanBootC-numStds(k)*stdC],'LineStyle','--','Color','k');hold on
        end
        plot(1:1:length(allDist{i,4}),allDist{i,4},'kx','MarkerSize',12);hold on
        plot(0,meanBootC,'k*');hold on
        constants=polyfit(1:1:length(allDist{i,4}),allDist{i,4},1);%fit best straight line to data points (least squares method)
        fitLine=constants(1)*[1:1:length(allDist{i,4})]+constants(2);
%         plot(0:1:numSessions-1,fitLine,'k:');hold on
        set(gca,'XLim',[0 numSessions]);
        set(gca,'XTick',0:5:numSessions);
        set(gca,'XTickLabel',0:5:numSessions);
        set(gca,'YLim',[0 1]);
        numTrials=grandTrialsList{rowInd,2}([grandTrialsList{rowInd,3}]==allDist{i,2});%number of trials for that session and that cell
        ptext=sprintf('%s  %s    # of trials per PSTH: %d',num2str(allDist{i,1}),num2str(allDist{i,2}),numTrials);
        text('Position',[0 1.05],'FontSize',9,'String',ptext);
        text('Position',[(numSessions)/1.2 0.17],'FontSize',9,'String',[char(8212),' within-session mean'],'Color','k');
        text('Position',[(numSessions)/1.2 0.13],'FontSize',9,'String',stdDevText,'Color','k');
        colInd=find(sessionNums==allDist{i,2});

        %calculate proportion of within-cell-between-session comparisons that fall within CIs:
        CItable(rowInd,colInd,1)=sum(allDist{i,4}(allDist{i,4}<=meanBootC+1.96*stdC)>=meanBootC-2*stdC)/length(allDist{i,4});
        %calculate proportion of between-cell-between-session comparisons that fall within CIs:
        CItable(rowInd,colInd,2)=sum(allDist{i,5}(allDist{i,5}<=meanBootC+1.96*stdC)>=meanBootC-2*stdC)/length(allDist{i,5});
        %calculate proportion of within-cell-between-session comparisons that fall within CIs:
        CItable(rowInd,colInd,3)=sum(allDist{i,4}(allDist{i,4}<=meanBootC+2.58*stdC)>=meanBootC-3*stdC)/length(allDist{i,4});
        %calculate proportion of between-cell-between-session comparisons that fall within CIs:
        CItable(rowInd,colInd,4)=sum(allDist{i,5}(allDist{i,5}<=meanBootC+2.58*stdC)>=meanBootC-3*stdC)/length(allDist{i,5});
        proportionText1=sprintf('wCellbSess  %.3f     (proportion of coef values within 95%% CI)',CItable(rowInd,colInd,1));
        proportionText2=sprintf('bCellbSess  %.3f',CItable(rowInd,colInd,2));
        proportionText3=sprintf('wCellbSess  %.3f     (proportion of coef values within 99%% CI)',CItable(rowInd,colInd,3));
        proportionText4=sprintf('bCellbSess  %.3f',CItable(rowInd,colInd,4));
        text('Position',[(numSessions)/2 1.05],'FontSize',9,'String',proportionText1,'Color','k');
        text('Position',[(numSessions)/2 1.01],'FontSize',9,'String',proportionText2,'Color','r');
        text('Position',[(numSessions)/2 0.97],'FontSize',9,'String',proportionText3,'Color','k');
        text('Position',[(numSessions)/2 0.93],'FontSize',9,'String',proportionText4,'Color','r');
            
        %distribution of coefs comparing current cell with OTHER cells and sessions:
        meanbCellC=mean(allDist{i,5});%between-cells
        stdbCellC=std(allDist{i,5});
        rowInd=find([grandTrialsList{:,1}]==allDist{i,1});
        all_sessions=grandTrialsList{rowInd,3};
        sessions=all_sessions(all_sessions~=allDist{i,2});
        all_sessions(all_sessions==3551)=355;all_sessions(all_sessions==3552)=355.5;
        sessions(sessions==3551)=355;sessions(sessions==3552)=355.5;
        line([0 numSessions],[meanbCellC meanbCellC],'LineStyle','-','Color','r');hold on
        for k=1:length(numStds)
            line([0 numSessions],[meanbCellC+numStds(k)*stdbCellC meanbCellC+numStds(k)*stdbCellC],'LineStyle','--','Color','r');hold on
            line([0 numSessions],[meanbCellC-numStds(k)*stdbCellC meanbCellC-numStds(k)*stdbCellC],'LineStyle','--','Color','r');hold on
        end
        text('Position',[(numSessions)/1.2 0.09],'FontSize',9,'String',[char(8212),' across-cells mean'],'Color','r');
        text('Position',[(numSessions)/1.2 0.05],'FontSize',9,'String',stdDevText,'Color','r');
        ch=[num2str(allDist{i,1}),'_',num2str(allDist{i,2})];
        imageName=[num2str(ch),'_',area,'_CI_PSTHs'];
        imageFolder=fullfile('F:','PL','xcorr',animal,subfolder);
        imagePath=fullfile(imageFolder,imageName);
        printtext=['print -dpng ',imagePath];
        set(gcf,'PaperPositionMode','auto')
        eval(printtext);
        allDist{i,7}=constants(1);%save the value of slope of best fit line to data points to .mat array
        close all
    end
end
matTexts=[{'skewedBoot'} {'CItable'} {'allDist'}];
for i=1:length(matTexts)
    matName=[matTexts{i},'_',area];
    matPath=fullfile('F:','PL','xcorr',animal,subfolder,matName);
    saveText=['save ',matPath,' ',matTexts{i}];
    eval(saveText);
end

