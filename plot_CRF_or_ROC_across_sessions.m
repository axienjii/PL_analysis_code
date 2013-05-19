function [slopeNeuro,c50,diffc50,minRate,maxRate,chSSE,yLimData,threshold82lower,threshold82higher]=plot_CRF_or_ROC_across_sessions(animal,area,analysisType,dataArray,chNum,numsessions,sessionSorted1,sampleContrast,testContrast,calculateTangent,plotDiffC50_30,excludeSessHighSSE,excludeOutliers,SSEMatPath,startEndTime,slC50MatPathname,slSigmaMultiple,c50SigmaMultiple,threshSigmaMultiple,rootFolder,plotLeastSquares)
if ischar(chNum)
    titleText=chNum;
    chNum=0;
else 
    titleText=num2str(chNum);
end
appendText=['_',num2str(sampleContrast)];
if ~strcmp(analysisType,'ROC_diff')
    fighandle1=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
    set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
end
chSSE=[];
slopeNeuro=[];
c50=[];
diffc50=[];
minRate=[];
maxRate=[];
yLimData=[];
threshold82lower=[];
threshold82higher=[];
for i=1:numsessions
    datavals=dataArray{i,3};%for channel of interest, for each session
    subplot(ceil(numsessions/5),5,i);
    if strcmp(analysisType,'CRF')
        xvals=testContrast(1):1:testContrast(end);
        plot(testContrast,datavals,'ok');
        hold on
        fittingFunction='n_r';
        options = optimset('Display','off','MaxFunEvals',10^6,'MaxIter',10^6,'TolFun',1.0E-6,'TolX',1.0E-6);
        X0=[max(datavals) 30 0.5 min(datavals)];
        if mean(datavals(1:3))>mean(datavals(12:14))%||chNum==13.2||chNum==24||chNum==42
            X0=[max(datavals) 30 -0.1 min(datavals)];%negative slope
        end
        [X]=fminsearch(fittingFunction,X0,options,testContrast,datavals);
        
        if calculateTangent==0
            slopeNeuro(1,i)=X(3);
        elseif calculateTangent==1
            slopeNeuro(1,i)=X(2)^X(3)*sampleContrast^(X(3)-1)/(30^X(3)+X(2)^X(3))^2;
        end
        c50(1,i)=X(2);
        minRate(1,i)=X(4);
        maxRate(1,i)=X(1);
        %     [X,fval]=fminsearch('weib_sim_min_max',X0,options,testContrast,datavals);
        %     slopeNeuro(1,i)=X(2);
        %     c50(1,i)=X(1).*(-log(0.5)).^(1/X(2));
        %     yvals=max(datavals)-(max(datavals)-min(datavals)).*exp(-((xvals./X(1
        %     )).^X(2)));
        if plotDiffC50_30==1
            diffc50(1,i)=abs(c50(1,i)-30);
        else
            diffc50=[];
        end
        yvals=X(1)*(xvals.^X(3)./(xvals.^X(3)+X(2).^X(3)))+X(4);
        plot(xvals,yvals,'Color','k');
        fitted_yvals=X(1)*(testContrast.^X(3)./(testContrast.^X(3)+X(2).^X(3)))+X(4);
        residuals=datavals-fitted_yvals;
        sseCRF=sum(residuals.^2);
        chSSE(i,1:2)=[i sseCRF];
    elseif strcmp(analysisType,'ROC')||strcmp(analysisType,'ROC_diff')||strcmp(analysisType,'ROC_diff2')
        xvals=testContrast(1):1:testContrast(end);
        if sum(datavals(1:3))<sum(datavals(end-2:end))||chNum==24&&sessionSorted1(1,i)==335||chNum==24&&sessionSorted1(1,i)==336||chNum==13
            X0=[2 30 0.2 0.1];
        elseif sum(datavals(1:3))>sum(datavals(end-2:end))
            X0=[-2 30 0.2 0.1];
        end
        options = optimset('Display','off','MaxFunEvals',10^4,'MaxIter',10^4,'TolFun',1.0E-6,'TolX',1.0E-6);
        X1=fminsearch(@fit_weibull,X0,options,testContrast,datavals,[],'least_square',[1 1 1 0],[10 100 1 0],[1 1 0 0],[-20 0 0 0]);
%         X=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 1],[2 100 1 0.2],[1 1 0 0],[-20 0 0 0]);
        X=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 0],[10 100 1 1],[1 1 0 0],[-10 0 0 0]);
        fitted_yvals=1-X(4)-X(3).*exp(-(testContrast./X(2)).^X(1));
        maxRate(1,i)=X(3);
        minRate(1,i)=X(4);
        residuals=datavals-fitted_yvals;
        sseCRF=sum(residuals.^2);
        chSSE(i,1:2)=[i sseCRF];
        if sseCRF>0.1%if the fit seems poor, try a variety of values for the upper and lower limits
            coefEsts2 = zeros(6,4);
            InitVar=1;
            for upperMax=[X1(3) X1(3)+0.1]
                for lowerMin=[X1(4) X1(4)+0.1]
                    [coefEsts2(InitVar,:)]=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 0],[20 100 upperMax 0],[1 1 0 1],[-20 0 0 lowerMin]);
                    fitted_yvals=1-coefEsts2(InitVar,4)-coefEsts2(InitVar,3).*exp(-(testContrast./coefEsts2(InitVar,2)).^coefEsts2(InitVar,1));
                    residuals=datavals-fitted_yvals;
                    sseCRFtemp(InitVar)=sum(residuals.^2);
                    InitVar=InitVar+1;
                end
            end
            [minSSE I]=min(sseCRFtemp);
            if minSSE<sseCRF
                X=coefEsts2(I,:);
                chSSE(i,1:2)=[i minSSE];
            end
        end
        if calculateTangent==0
            slopeNeuro(1,i)=X(1);
        elseif calculateTangent==1
            slopeNeuro(1,i)=X(1)*X(3)*exp(-(sampleContrast/X(2))^X(1))*sampleContrast^(X(1)-1)*(1/X(2))^X(1);
        end
        c50(1,i)=real(X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1)));
        hold on
        if plotDiffC50_30==1
            diffc50(1,i)=abs(c50(1,i)-30);
        end
        if strcmp(analysisType,'ROC_diff2')
            line(c50(1,i),0:0.01:1,'Color','b');
        else
            line(c50(1,i),0:0.01:1,'Color','r');
        end
        yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
%         yvals=(yvals-.005)/.99;
%         yvals=log(yvals);
        if strcmp(analysisType,'ROC_diff2')
            plot(xvals,yvals,'Color','b');
        elseif strcmp(analysisType,'ROC_diff')
            plot(xvals,yvals,'Color','r');
        else
            if plotLeastSquares==1
                leastSquares_yvals=1-X1(4)-X1(3).*exp(-(xvals./X1(2)).^X1(1));
                plot(xvals,leastSquares_yvals,'Color','b');%calculated using least squares minimization instead of maximum likelihood estimation
            end
            plot(xvals,yvals,'Color','k');
        end
        hold on
        if strcmp(analysisType,'ROC_diff2')
%             plot(testContrast,datavals,'ob');%empty markers
            plot(testContrast,datavals,'ob','MarkerFaceColor','b');
        elseif strcmp(analysisType,'ROC_diff')
            %read R-values for each session and condition
            matName=[num2str(chNum),'_Rp_ROCdiff_',area];
            matPathName=fullfile(rootFolder,'PL','ROC',animal,'ROC_diff',matName);
            loadText=['load ',matPathName,' p r'];
            eval(loadText);
            yLimVals=get(gca,'ylim');
            xLimVals=get(gca,'xlim');
            unitSpace=(yLimVals(2)-yLimVals(1))/30;
            colmapText=colormap(jet(size(testContrast,2)));
            colmapText=[colmapText(1,:);132/255 22/255 216/255;202/255 65/255 223/255;colmapText(3:7,:);157/255 212/255 61/255;colmapText(10:12,:);178/255 111/255 12/255;colmapText(end,:)];
            for cond=1:length(testContrast)
%                 plot(testContrast(cond),datavals(cond),'Color','r','Marker','o','LineStyle','none');%empty markers
                plot(testContrast(cond),datavals(cond),'Color','r','Marker','o','LineStyle','none','MarkerFaceColor','r');
%                 plot(testContrast(cond),(r(i,cond)+1)/2,'Color',colmapText(cond,:),'Marker','o','LineStyle','none');%different colours for each cond when plotting R values
                plot(testContrast(cond),(r(i,cond)+1)/2,'Color','g','Marker','o','LineStyle','none','MarkerFaceColor','g');
%                 [AX,H1,H2] = plotyy(testContrast(cond),datavals(cond),testContrast(cond),(r(i,cond)+1)/2,'plot');
%                 set(H1,'Color','r','Marker','o','LineStyle','none');
%                 set(H3,'Color',colmapText(cond,:),'Marker','o','LineStyle','none');
%                 set(AX1,'XTickLabel',[]);
%                 set(AX1,'YTickLabel',[]);
%                 set(AX2,'XTickLabel',[]);
%                 set(AX2,'YTickLabel',[]);
                text('Position',[xLimVals(2)+(xLimVals(2)-xLimVals(1))/25 yLimVals(1)+unitSpace*cond*2.2],'FontSize',9,'String',sprintf('%.2f',r(i,cond)),'Color',colmapText(cond,:));
            end
            ylim([0 1]);
%             set(AX1,'ycolor','k')
%             set(AX2,'ycolor','k')
%             set(AX1,'YLim',[0 1]);
%             set(AX1,'XTickLabel',[0 20 40 60]);
%             set(AX1,'YTick',[0 0.5 1],'YTickLabel',[0 0.5 1]);
%             set(AX2,'YLim',[0 1]);%R vals have been transformed so range is no longer [-1 1] but [0 1]
%             set(AX2,'YTick',[0 0.5 1],'YTickLabel',[-1 0 1]);%relabel y-axis ticks
        else
            plot(testContrast,datavals,'ok');
        end
    elseif strcmp(analysisType,'NVP')||strcmp(analysisType,'psycho')||strcmp(analysisType,'psycho_param')
        colHighLow='rb';
        maxDiff=[sampleContrast 100-sampleContrast];%maximum possible difference between sample contrast and extremes of contrast range (i.e. [0 100])
        for higherOrLower=1:2%separate curve fits for data where test contrast is below or above sample contrast
            ind=find(testContrast<sampleContrast);
            ind=ind(end);
            if higherOrLower==1%lower contrast than sample
                testContrastSplit(higherOrLower)={[0 fliplr(sampleContrast-testContrast(1:ind))]};%calculate absolute difference between sample and test contrasts
                datavalsSplit(higherOrLower)={[0.5 fliplr(1-datavals(1:ind))]};
            elseif higherOrLower==2%higher contrast than sample 
                testContrastSplit(higherOrLower)={[0 testContrast(ind+1:end)-sampleContrast]};  
                datavalsSplit(higherOrLower)={[0.5 datavals(ind+1:end)]};    
            end
            xvals=0:0.1:60;%hard coded 0:0.1:testContrastSplit{higherOrLower}(end)
            plot([testContrastSplit{higherOrLower}],[100*datavalsSplit{higherOrLower}],'Marker','o','LineStyle','none','Color',colHighLow(higherOrLower));
            hold on
            options = optimset('Display','off','MaxFunEvals',10^6,'MaxIter',10^6,'TolFun',1.0E-6,'TolX',1.0E-6);
            if strcmp(analysisType,'psycho')||strcmp(analysisType,'psycho_param')
                if strcmp(analysisType,'psycho_param')
                    maxSlope=1.8;
                    stepSize=0.2;
                else
                    maxSlope=9.5;
                    stepSize=2;
                end
                coefEsts2 = zeros(20,2);
                fval2 = zeros(20,1);
                InitVar = 0;                
                for iC50 = 1 : 4                    
                    nC50 = iC50 * 5;                    
                    for iSlope = 0.5 : stepSize : maxSlope                        
                        Slope = iSlope;
                        InitVar = InitVar + 1;
                        X0=[nC50 Slope];
                        if strcmp(analysisType,'psycho')
                            [coefEsts2(InitVar,:),fval2(InitVar)]=fminsearch(@weibull_pointfive_one,X0,options,testContrastSplit{higherOrLower},datavalsSplit{higherOrLower},[1 1],[maxDiff(higherOrLower) 10],[1 1],[0 0],'mle');%place upper limit on threshold value, depending on whether examining higher or lower test contrast conds
                        elseif strcmp(analysisType,'psycho_param')
                            [coefEsts2(InitVar,:),fval2(InitVar)]=fminsearch(@weibull_behav_param,X0,options,testContrastSplit{higherOrLower},datavalsSplit{higherOrLower},[1 1],[maxDiff(higherOrLower) 10],[1 1],[0 0],1-datavalsSplit{higherOrLower}(end),'mle');%place upper limit on threshold value, depending on whether examining higher or lower test contrast conds
                        end
                    end
                end
                [dummy I]=min(fval2);
                coefEsts=coefEsts2(I,:);
                X=coefEsts;
%                 if strcmp(analysisType,'psycho')
%                     X=fminsearch(@weibull_pointfive_one,coefEsts,options,testContrastSplit{higherOrLower},datavalsSplit{higherOrLower},[1 1],[maxDiff(higherOrLower) 10],[1 1],[0 0]);%place upper limit on threshold value, depending on whether examining higher or lower test contrast conds
%                 elseif strcmp(analysisType,'psycho_param')
%                     X=fminsearch(@weibull_behav_param,X0,options,testContrastSplit{higherOrLower},datavalsSplit{higherOrLower},[1 1],[maxDiff(higherOrLower) 1.8],[1 1],[0 0],1-datavalsSplit{higherOrLower}(end));%includes proportion of erroneous responses for condition with highest test contrast
%                 end
            else
                if sum(datavalsSplit{higherOrLower}(1:3))<sum(datavalsSplit{higherOrLower}(end-2:end))||chNum==24&&sessionSorted1(1,i)==335||chNum==24&&sessionSorted1(1,i)==336||chNum==13
                    X0=[10 1];
                elseif sum(datavalsSplit{higherOrLower}(1:3))>sum(datavalsSplit{higherOrLower}(end-2:end))
                    X0=[10 -1];
                end
                X=fminsearch(@weibull_pointfive_one,X0,options,testContrastSplit{higherOrLower},datavalsSplit{higherOrLower},[1 1],[maxDiff(higherOrLower) 10],[1 1],[0 0],'mle');%place upper limit on threshold value, depending on whether examining higher or lower test contrast conds
            end
            hold on
            if strcmp(analysisType,'psycho')
            if higherOrLower==1%lower contrast than sample
                threshold82lower(1,i)=X(1);%value of 81.161 is obtained when the contrast (x) is equal to the threshold (1-0.5*exp(-1) = 0.8161). value of 0.36 is from 2*(1-0.82)
            elseif higherOrLower==2%higher contrast than sample
                threshold82higher(1,i)=X(1);
            end
            elseif strcmp(analysisType,'psycho_param')
                if higherOrLower==1%lower contrast than sample
                    threshold82lower(1,i)=X(1).*(-log(1-(1-0.5*exp(-1)-0.5)./(0.5-(1-datavalsSplit{higherOrLower}(end))))).^(1/X(2));%value of 81.161 is obtained when the contrast (x) is equal to the threshold (1-0.5*exp(-1) = 0.8161). value of 0.36 is from 2*(1-0.82)
                    lineTypeLower(1,i)={'--'};
                    if ~isreal(threshold82lower(1,i))||threshold82lower(1,i)>sampleContrast
                        threshold82lower(1,i)=sampleContrast;
                        lineTypeLower(1,i)={'-.'};
                    end
                elseif higherOrLower==2%higher contrast than sample
                    threshold82higher(1,i)=X(1).*(-log(1-(1-0.5*exp(-1)-0.5)./(0.5-(1-datavalsSplit{higherOrLower}(end))))).^(1/X(2));
                    lineTypeHigher(1,i)={'--'};
                    if ~isreal(threshold82higher(1,i))||threshold82higher(1,i)>100-sampleContrast
                        threshold82higher(1,i)=100-sampleContrast;
                        lineTypeHigher(1,i)={'-.'};
                    end
                end
            end
            if strcmp(analysisType,'psycho_param')
                yvals=100*(0.5+(0.5-(1-datavalsSplit{higherOrLower}(end))).*(1-exp(-(xvals./X(1)).^X(2))));
            else
                yvals=100*(1-0.5.*exp(-(xvals./X(1)).^X(2)));
            end
            plot(xvals,yvals,'Color',colHighLow(higherOrLower));  
            if strcmp(analysisType,'psycho_param')
                fitted_yvals=100*(0.5+(0.5-(1-datavalsSplit{higherOrLower}(end))).*(1-exp(-(testContrastSplit{higherOrLower}./X(1)).^X(2))));
            else
                fitted_yvals=100*(1-0.5.*exp(-(testContrastSplit{higherOrLower}./X(1)).^X(2)));
            end
            residuals=datavalsSplit{higherOrLower}-fitted_yvals/100;
            sseCRF=sum(residuals.^2);
            chSSE(i,1:2)=[i sseCRF];
        end
    end
    if strcmp(analysisType,'ROC')||strcmp(analysisType,'CRF')
        ylim([0,max(datavals)]);
%         ylim([0 1]);
    end
    xlim([0 max(testContrast)]);
    subplottitle=num2str(i);
    title(subplottitle);
    [yLimData(i,1:2)]=get(gca,'YLim');
end
if strcmp(analysisType,'NVP')||strcmp(analysisType,'psycho')||strcmp(analysisType,'psycho_param')
    for i=1:numsessions
        subplot(ceil(numsessions/5),5,i);
        for higherOrLower=1:2
            if strcmp(analysisType,'psycho_param')
                if higherOrLower==1%lower contrast than sample
                    set(gca,'YLim',yLimData(i,1:2));
                    plot([threshold82lower(1,i) threshold82lower(1,i)],[yLimData(i,1) yLimData(i,2)],'Color',colHighLow(higherOrLower),'LineStyle',lineTypeLower{1,i});
                elseif higherOrLower==2%higher contrast than sample
                    plot([threshold82higher(1,i) threshold82higher(1,i)],[yLimData(i,1) yLimData(i,2)],'Color',colHighLow(higherOrLower),'LineStyle',lineTypeHigher{1,i});
                end
            else
                if higherOrLower==1%lower contrast than sample
                    set(gca,'YLim',yLimData(i,1:2));
                    plot(threshold82lower(1,i),yLimData(i,1):(yLimData(i,2)-yLimData(i,1))/100:yLimData(i,2),'Color',colHighLow(higherOrLower));
                elseif higherOrLower==2%higher contrast than sample
                    plot(threshold82higher(1,i),yLimData(i,1):(yLimData(i,2)-yLimData(i,1))/100:yLimData(i,2),'Color',colHighLow(higherOrLower));
                end
            end
        end
        xlim([0 max([testContrastSplit{1,1} testContrastSplit{1,2}])+10]);
    end
elseif strcmp(analysisType,'ROC')||strcmp(analysisType,'CRF')||strcmp(analysisType,'ROC_diff2')
    set_ylim_across_sessions(titleText,numsessions,sampleContrast,testContrast,c50,chSSE,yLimData,analysisType);
end
imageTitleText=titleText;
if find(titleText==' ')
    imageTitleText(imageTitleText==' ')='_';
end    
if strcmp(analysisType,'ROC_diff')
    imagename=[imageTitleText,appendText,startEndTime,'_',area,'_dualPlots'];
end
if excludeSessHighSSE==0&&~strcmp(analysisType,'ROC_diff')&&~strcmp(analysisType,'ROC_diff2')
    saveText=['save ',SSEMatPath,' chSSE'];
    eval(saveText);
    imagename=[imageTitleText,appendText,startEndTime,'_',area];
    if strcmp(analysisType,'psycho_param')%include all sessions for behavioural paper
        sessionSorted2=sessionSorted1;
        saveText=['save ',slC50MatPathname,' sessionSorted2 threshold82lower threshold82higher'];
        eval(saveText)
    end
elseif excludeSessHighSSE==1
    if excludeOutliers==0
        imagename=[imageTitleText,appendText,startEndTime,'_',area,'_goodSSE'];
        if strcmp(analysisType,'NVP')
            saveText=['save ',slC50MatPathname,'.mat sessionSorted1 threshold82lower threshold82higher'];
            eval(saveText)
        elseif strcmp(analysisType,'psycho')||strcmp(analysisType,'psycho_param')
            sessionSorted2=sessionSorted1;
            saveText=['save ',slC50MatPathname,' sessionSorted2 threshold82lower threshold82higher'];
            eval(saveText)
        elseif strcmp(analysisType,'ROC')||strcmp(analysisType,'CRF')
            saveText=['save ',slC50MatPathname,'.mat sessionSorted1 slopeNeuro c50'];
            eval(saveText)
        end
    elseif excludeOutliers==1
        if strcmp(analysisType,'ROC')||strcmp(analysisType,'CRF')
            imagename=[imageTitleText,appendText,startEndTime,'_',area,'_goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
        elseif strcmp(analysisType,'NVP')||strcmp(analysisType,'psycho')
            imagename=[imageTitleText,appendText,startEndTime,'_',area,'_goodSSE_no_outliers_th',num2str(threshSigmaMultiple)];
            sessionSorted2=sessionSorted1;
            if strcmp(analysisType,'NVP')
                saveText=['save ',slC50MatPathname,' sessionSorted1 threshold82lower threshold82higher'];
            elseif strcmp(analysisType,'psycho')||strcmp(analysisType,'psycho_param')
                saveText=['save ',slC50MatPathname,' sessionSorted2 threshold82lower threshold82higher'];
            end
            eval(saveText)
        end
    end
end
subFolder=['neurometric_',analysisType];
if ~strcmp(analysisType,'ROC_diff2')
    pathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder,imagename);
    if strcmp(analysisType,'ROC_diff')
        pathname=fullfile(rootFolder,'PL','ROC',animal,'ROC_diff',imagename);
    end
    folder=fullfile(rootFolder,'PL',analysisType,animal,subFolder);
    if ~exist(folder,'dir')
        mkdir(folder)
    end
    printtext=sprintf('print -dpng %s.png',pathname);
    set(gcf,'PaperPositionMode','auto')
    eval(printtext);
end
% appendData(1:numsessions,1)=chNum;
% appendData(1:numsessions,2)=sessionSorted1;
% appendCRF(1:numsessions,3:length(testContrast)+2)=VALUES;
% appendData=[appendData,threshold82lower',threshold82lower'];