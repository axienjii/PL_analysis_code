function [slopeNeuro,c50,diffc50,minRate,maxRate,chSSE,yLimData]=plot_CRF_or_ROC_across_sessions(animal,area,analysisType,dataArray,chNum,numsessions,sessionSorted1,sampleContrast,testContrast,calculateTangent,plotDiffC50_30,excludeSessHighSSE,excludeOutliers,SSEMatPath,startEndTime,slC50MatPathname,slSigmaMultiple,c50SigmaMultiple)
if ischar(chNum)
    titleText=chNum;
    chNum=0;
else 
    titleText=num2str(chNum);
end
appendText=['_',num2str(sampleContrast)];
fighandle1=  figure('Color',[1,1,1],'Units','Normalized','Position',[0.1, 0.1, 0.8, 0.8]); %
set(fighandle1, 'PaperUnits', 'centimeters', 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPosition', [0.63452 0.63452 6.65 3.305]);
for i=1:numsessions
    datavals=dataArray{i,3};%for channel of interest, for each session
    xvals=testContrast(1):1:testContrast(end);
    subplot(ceil(numsessions/5),5,i);
    plot(testContrast,datavals,'ok');
    hold on
    if strcmp(analysisType,'CRF')
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
            slopeNeuro(1,i)=X(2)^X(3)*30^(X(3)-1)/(30^X(3)+X(2)^X(3))^2;
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
    elseif strcmp(analysisType,'ROC')
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
            slopeNeuro(1,i)=X(1)*X(3)*exp(-(30/X(2))^X(1) )*30^(X(1)-1)*(1/X(2))^X(1);
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
    end
    residuals=datavals-fitted_yvals;
    sseCRF=sum(residuals.^2);
    chSSE(i,1:2)=[i sseCRF];
    plot(xvals,yvals,'r');
    %         ylim([0,max(datavals)]);
    % ylim([0 1]);
    xlim([0 max(testContrast)]);
    subplottitle=num2str(i);
    title(subplottitle);
    [yLimData(i,1:2)]=get(gca,'YLim');
end
set_ylim_across_sessions(titleText,numsessions,sampleContrast,testContrast,c50,chSSE,yLimData,analysisType);
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
        saveText=['save ',slC50MatPathname,'.mat sessionSorted1 slopeNeuro c50'];
        eval(saveText)
    elseif excludeOutliers==1
        imagename=[imageTitleText,appendText,startEndTime,'_',area,'_goodSSE_no_outliers_sl',num2str(slSigmaMultiple),'_C50',num2str(c50SigmaMultiple)];
    end
end
subFolder=['neurometric_',analysisType];
pathname=fullfile('F:','PL',analysisType,animal,subFolder,imagename);
folder=fullfile('F:','PL',analysisType,animal,subFolder);
if ~exist(folder,'dir')
    mkdir(folder)
end
printtext=sprintf('print -dpng %s.png',pathname);
eval(printtext);
appendData(1:numsessions,1)=chNum;
appendData(1:numsessions,2)=sessionSorted1;
% appendCRF(1:numsessions,3:length(testContrast)+2)=VALUES;
if plotDiffC50_30==1
    appendData=[appendData,slopeNeuro',c50',diffc50',minRate',maxRate'];
else
    appendData=[appendData slopeNeuro' c50' minRate' maxRate'];
end