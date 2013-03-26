function read_bj_SE_xcorr5skewslope(animal,area,channels)
%Written by Xing 18/10/10
%Reads values of slope and intercept of fitted data (consisting of correlation coefficients between the session
%and other sessions for that cell) with a straight line, that were
%calcualted in read_blanco_SE_xcorr5skew, using the polyfit 
%function. Then, examines values of
%slopes for each cell across sessions to check whether the slope goes from
%more negative to more positive, over time.
%
%When you examine the degree to which PSTH data from each session
%correlates to that from all other sessions, you would expect to find that PSTH values from
%each session correlate best with data from sessions close to it in time-
%for example, data from an early session would likely be most highly-correlated with data from other
%early sessions, whereas data from a late session would likely be most highly-correlated
%with that from other late sessions. Thus, the pattern tends to be that the
%degree of correlation between data from a session in the middle of the
%training period (e.g. session 327 or 328) with data from the rest of the
%sessions, is higher than that between data from an early or late session
%with all other sessions.          
%This is shown by the slope of the best-fit line to data points marked by
%crosses in the following figures. If the value of the slope is negative for
%early sessions, then this indicates that levels of correlation are high between
%the early session with other early sessions, then get progressively lower.
%If the value of the slope becomes less negative, and reaches 0 around the
%middle of the training period, this indicates that each middle session is as
%well correlated to early as to late sessions. If the value of the slope then
%becomes progressively positive for later sessions, this indicates that the
%correlation between a late session with early sessions is low, while its
%correlation with other late sessions is high.   
%Saves image of change in slope for each channel, naming conventions:
%-xcorrslope_min8datapoints: at least 8 data points per plot for that
%channel to be included
%-xcorrslope_min8datap_redline: additionally, draws red line calculated during linear
%regression fitting, for each channel.

%compile values of slopes across sessions for each channel and check
%whether, as hypothesised, the slopes indeed go from more negative to more
%positive, for each channel:

minusSpontan=0;
sigma=8;
%combine PSTH activity values across channels and sessions:
if nargin<3 || isempty(channels)
    channels = main_channels(animal,area);
end
if minusSpontan==1
    subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_mspontan_',area];%folder for stimulus-evoked responses minus spontaneous activity levels
elseif minusSpontan==0
    subfolder=['PSTH45_images_sm',num2str(sigma*2),'ms_wspontan_',area];%folder for stimulus-evoked responses without any subtraction of spontaneous activity levels
end
matName=['allDist_',area];
matPath=fullfile('F:','PL','xcorr',animal,subfolder,matName);
loadText=['load ',matPath,' allDist'];
eval(loadText);

allSlopes=[];
for channelNum=1:length(channels)
    slopes=[];
    slope_sessions=[];
    for row=1:size(allDist,1)
        if ~isempty(allDist{row,4})
            if allDist{row,1}==channels(channelNum);
                slopes=[slopes allDist{row,7}];
                session=allDist{row,2};
                if session==3551
                    session=355;
                end
                if session==3552
                    session=355.5;
                end
                slope_sessions=[slope_sessions session];
            end
        end
    end
    if ~isempty(slopes)
        if length(slopes)>5%adjust this to vary criterion for inclusion of channel based on # of data points available
            allSlopes=[allSlopes;channels(channelNum) {slopes} {slope_sessions}];
        end
    end
end

fighandle1=figure('Color',[1,1,1],'Units','Normalized','Position',[0.02, 0.02, 0.8, 0.9]);
set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
orient landscape
for row=1:size(allSlopes,1)
    xSessions=1:length(allSlopes{row,2});
    subplot(ceil(size(allSlopes,1)/5),5,row);
    plot(xSessions,allSlopes{row,2},'kx');hold on
    constants=polyfit(xSessions,allSlopes{row,2},1);%fit best straight line to data points (least squares method)
    fitLine=constants(1)*xSessions+constants(2);
    plot(xSessions,fitLine,'r:');hold on
    a=[xSessions' allSlopes{row,2}'];%check for change in slope with time- expect a +ve correlation, if any
    [coefficients1 p1]=corrcoef(a);   
    allSlopes(row,4)={coefficients1(1,2)};
    allSlopes(row,5)={p1(1,2)};
    subplottitle=num2str(allSlopes{row,1});
    if p1(1,2)<0.05
        subplottitle=[subplottitle,'*'];
    end
    title(subplottitle);
    set(gca,'XLim',[0 length(allSlopes{row,2})]);
    set(gca,'XTick',[0 length(allSlopes{row,2})]);
    set(gca,'XTickLabel',[0 length(allSlopes{row,2})]);
    set(gca,'YLim',[min(allSlopes{row,2}) max(allSlopes{row,2})]);
    if ceil(min(allSlopes{row,2})*1000)/1000<floor(max(allSlopes{row,2})*1000)/1000
        set(gca,'YTick',[ceil(min(allSlopes{row,2})*1000)/1000 floor(max(allSlopes{row,2})*1000)/1000]);
        set(gca,'YTickLabel',[ceil(min(allSlopes{row,2})*1000)/1000 floor(max(allSlopes{row,2})*1000)/1000]);
    else
        set(gca,'YTick',[min(allSlopes{row,2}*1000)/1000 max(allSlopes{row,2})*1000/1000]);
        set(gca,'YTickLabel',[min(allSlopes{row,2})*1000/1000 max(allSlopes{row,2}*1000)/1000]);
    end
    ptext1=[' R=',num2str(coefficients1(1,2))];
    text('Position',[0 max(allSlopes{row,2})-0.1*(max(allSlopes{row,2})-min(allSlopes{row,2}))],'FontSize',7,'String',ptext1);
    ptext2=[' p=',num2str(p1(1,2))];
    text('Position',[0 max(allSlopes{row,2})-0.2*(max(allSlopes{row,2})-min(allSlopes{row,2}))],'FontSize',7,'String',ptext2);
    if row==1
        ptext=sprintf('%s    %s',matPath,area);
        orient landscape
        yLimVals=get(gca,'YLim');
        xLimVals=get(gca,'XLim');
        text('Position',[xLimVals(1)-0.1*(xLimVals(2)-xLimVals(1)) yLimVals(2)+0.3*(yLimVals(2)-yLimVals(1))],'FontSize',10,'String',ptext);
    end
end
imageName=['xcorrslope_',area,'_PSTHc'];
imageFolder=fullfile('F:','PL','xcorr',animal,subfolder);
imagePath=fullfile(imageFolder,imageName);
printtext=['print -dpng ',imagePath];
set(gcf,'PaperPositionMode','auto')
eval(printtext);
close all