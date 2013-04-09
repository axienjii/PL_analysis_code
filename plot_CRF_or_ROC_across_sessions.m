function [slopeNeuro,c50,diffc50,minRate,maxRate,chSSE,yLimData,threshold82lower,threshold82higher]=plot_CRF_or_ROC_across_sessions(animal,area,analysisType,dataArray,chNum,numsessions,sessionSorted1,sampleContrast,testContrast,calculateTangent,plotDiffC50_30,excludeSessHighSSE,excludeOutliers,SSEMatPath,startEndTime,slC50MatPathname,slSigmaMultiple,c50SigmaMultiple,threshSigmaMultiple,rootFolder)
if ischar(chNum)
    titleText=chNum;
    chNum=0;
else 
    titleText=num2str(chNum);
end
appendText=['_',num2str(sampleContrast)];
fighandle1=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
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
        fitted_yvals=X(1)*(testContrast.^X(3)./(testContrast.^X(3)+X(2).^X(3)))+X(4);
        residuals=datavals-fitted_yvals;
        sseCRF=sum(residuals.^2);
        chSSE(i,1:2)=[i sseCRF];
    elseif strcmp(analysisType,'ROC')
        xvals=testContrast(1):1:testContrast(end);
        plot(testContrast,datavals,'ok');
        hold on
        if sum(datavals(1:3))<sum(datavals(end-2:end))||chNum==24&&sessionSorted1(1,i)==335||chNum==24&&sessionSorted1(1,i)==336||chNum==13
            X0=[2 30 0.2 0.1];
        elseif sum(datavals(1:3))>sum(datavals(end-2:end))
            X0=[-2 30 0.2 0.1];
        end
        options = optimset('Display','off','MaxFunEvals',10^6,'MaxIter',10^6,'TolFun',1.0E-6,'TolX',1.0E-6);
        X=fminsearch(@fit_weibull,X0,options,testContrast,datavals,[],'least_square',[1 1 0 0],[20 100 0 0],[1 1 0 0],[-20 0 0 0]);
        if calculateTangent==0
            slopeNeuro(1,i)=X(1);
        elseif calculateTangent==1
            slopeNeuro(1,i)=X(1)*X(3)*exp(-(sampleContrast/X(2))^X(1) )*sampleContrast^(X(1)-1)*(1/X(2))^X(1);
        end
        c50(1,i)=real(X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1)));
        plot(testContrast,datavals,'ok');
        hold on
        if plotDiffC50_30==1
            diffc50(1,i)=abs(c50(1,i)-30);
        end
        line(c50(1,i),0:0.01:1,'Color','r');
        yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
        fitted_yvals=1-X(4)-X(3).*exp(-(testContrast./X(2)).^X(1));
        minRate(1,i)=X(3);
        maxRate(1,i)=X(4);
        residuals=datavals-fitted_yvals;
        sseCRF=sum(residuals.^2);
        chSSE(i,1:2)=[i sseCRF];
    elseif strcmp(analysisType,'NVP')||strcmp(analysisType,'psycho')
        colHighLow='rb';
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
            if sum(datavalsSplit{higherOrLower}(1:3))<sum(datavalsSplit{higherOrLower}(end-2:end))||chNum==24&&sessionSorted1(1,i)==335||chNum==24&&sessionSorted1(1,i)==336||chNum==13
                X0=[10 1];
            elseif sum(datavalsSplit{higherOrLower}(1:3))>sum(datavalsSplit{higherOrLower}(end-2:end))
                X0=[10 -1];
            end
            X=fminsearch(@weibull_pointfive_one,X0,options,testContrastSplit{higherOrLower},datavalsSplit{higherOrLower},[1 1],[60 5],[1 1],[0 0]);
            hold on
            if higherOrLower==1%lower contrast than sample
                threshold82lower(1,i)=X(1);%value of 81.161 is obtained when the contrast (x) is equal to the threshold (1-0.5*exp(-1) = 0.8161). value of 0.36 is from 2*(1-0.82)
            elseif higherOrLower==2%higher contrast than sample
                threshold82higher(1,i)=X(1);
            end
            yvals=100*(1-0.5.*exp(-(xvals./X(1)).^X(2)));
            plot(xvals,yvals,'Color',colHighLow(higherOrLower));         
            fitted_yvals=100*(1-0.5.*exp(-(testContrastSplit{higherOrLower}./X(1)).^X(2)));
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
if strcmp(analysisType,'NVP')||strcmp(analysisType,'psycho')
    for i=1:numsessions
        subplot(ceil(numsessions/5),5,i);
        for higherOrLower=1:2
            if higherOrLower==1%lower contrast than sample
                set(gca,'YLim',yLimData(i,1:2));
                plot(threshold82lower(1,i),yLimData(i,1):(yLimData(i,2)-yLimData(i,1))/100:yLimData(i,2),'Color',colHighLow(higherOrLower));
            elseif higherOrLower==2%higher contrast than sample
                plot(threshold82higher(1,i),yLimData(i,1):(yLimData(i,2)-yLimData(i,1))/100:yLimData(i,2),'Color',colHighLow(higherOrLower));
            end
        end
    end
elseif strcmp(analysisType,'ROC')||strcmp(analysisType,'CRF')
    set_ylim_across_sessions(titleText,numsessions,sampleContrast,testContrast,c50,chSSE,yLimData,analysisType);
end
imageTitleText=titleText;
if find(titleText==' ')
    imageTitleText(imageTitleText==' ')='_';
end    
if excludeSessHighSSE==0
    saveText=['save ',SSEMatPath,' chSSE'];
    eval(saveText);
    imagename=[imageTitleText,appendText,startEndTime,'_',area];
elseif excludeSessHighSSE==1
    if excludeOutliers==0
        imagename=[imageTitleText,appendText,startEndTime,'_',area,'_goodSSE'];
        if strcmp(analysisType,'NVP')
            saveText=['save ',slC50MatPathname,'.mat sessionSorted1 threshold82lower threshold82higher'];
            eval(saveText)
        elseif strcmp(analysisType,'psycho')
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
            elseif strcmp(analysisType,'psycho')
                saveText=['save ',slC50MatPathname,' sessionSorted2 threshold82lower threshold82higher'];
            end
            eval(saveText)
        end
    end
end
subFolder=['neurometric_',analysisType];
pathname=fullfile(rootFolder,'PL',analysisType,animal,subFolder,imagename);
folder=fullfile(rootFolder,'PL',analysisType,animal,subFolder);
if ~exist(folder,'dir')
    mkdir(folder)
end
printtext=sprintf('print -dpng %s.png',pathname);
set(gcf,'PaperPositionMode','auto')
eval(printtext);
% appendData(1:numsessions,1)=chNum;
% appendData(1:numsessions,2)=sessionSorted1;
% appendCRF(1:numsessions,3:length(testContrast)+2)=VALUES;
% appendData=[appendData,threshold82lower',threshold82lower'];