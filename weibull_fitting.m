function [slopeNeuro,c50,diffc50,minRate,maxRate,chSSE,threshold82higher,threshold82lower]=weibull_fitting(datavals,sampleContrast,testContrast,ROCanalysisType,i,slopeNeuro,chSSE,c50,minRate,maxRate,diffc50,plotDiffC50_30,calculateTangent,useISI,threshold82higher,threshold82lower)
alex_fit=0;
xvals=testContrast(1):1:testContrast(end);
if useISI==0
    if sum(datavals(1:3))<=sum(datavals(end-2:end))
        X0=[5 30 0.5 0.1];
    elseif sum(datavals(1:3))>sum(datavals(end-2:end))
        X0=[-3 30 0.5 0.1];
    end
elseif useISI==1
    if sum(datavals(1:3))<=sum(datavals(end-2:end))
        X0=[5 30 0.5 0];
    elseif sum(datavals(1:3))>sum(datavals(end-2:end))
        X0=[-3 30 0.5 0.5];
    end
end
options = optimset('Display','off','MaxFunEvals',10^4,'MaxIter',10^4,'TolFun',1.0E-6,'TolX',1.0E-6);
if alex_fit==1
    X0=[30 1];
    [X,fval(1)]=fminsearch(@weibull_zero_one,[X0],[],testContrast,datavals);
    fitted_yvals=1-exp(-((testContrast/X(1)).^X(2)));
    maxRate(1,i)=max(fitted_yvals);
    minRate(1,i)=min(fitted_yvals);
    residuals=datavals-fitted_yvals;
    sseCRF=sum(residuals.^2);
    chSSE(i,1:2)=[i sseCRF];
    if calculateTangent==0
        slopeNeuro(1,i)=X(2);
    elseif calculateTangent==1
        slopeNeuro(1,i)=X(2)*exp(-(sampleContrast/X(1))^X(2))*sampleContrast^(X(2)-1)*(1/X(1))^X(2);
    end
    
    yvals=1-exp(-((xvals/X(1)).^X(2)));
    xvalsFine=testContrast(1):0.01:testContrast(end);
    yvalsFine=1-exp(-((xvalsFine/X(1)).^X(2)));
    if useISI==0
        if max(yvals)<0.5%out of range
            c50(1,i)=100;
        elseif min(yvals)>0.5
            c50(1,i)=0;
        else
            diffTemp=yvalsFine-0.5;
            [tempVal columnInd]=min(abs(diffTemp));
            c50(1,i)=xvalsFine(columnInd);
            %                 c50(1,i)=real(X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1)));
        end
        threshold82higher=[];
    elseif useISI==1%read C50/PNE as being point where AUROC value is 0.75, instead of 0.5
        if max(yvals)<0.75%out of range
            c50(1,i)=100;
        elseif min(yvals)>0.75
            c50(1,i)=0;
        else
            diffTemp=yvalsFine-0.75;
            [tempVal columnInd]=min(abs(diffTemp));
            c50(1,i)=xvalsFine(columnInd);
            %                 c50(1,i)=real(X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1)));
        end
        diffTemp=[];
        if slopeNeuro(1,i)>=0&&max(yvals)>=0.82%stimulus-evoked excitation
            diffTemp=yvalsFine-0.82;%single neurometric threshold
        elseif slopeNeuro(1,i)<0&&min(yvals)<=1-0.82%stimulus-evoked suppression
            diffTemp=yvalsFine-(1-0.82);%single neurometric threshold
        end
        if ~isempty(diffTemp)
            if max(yvals)>=0.82&&min(yvals)<=0.82||max(yvals)>=1-0.82&&min(yvals)<=1-0.82
                [tempVal columnInd]=min(abs(diffTemp));
                threshold82higher(1,i)=xvalsFine(columnInd);
            else
                threshold82higher(1,i)=NaN;
            end
        else
            threshold82higher(1,i)=NaN;
        end
    end
else
    X1=fminsearch(@fit_weibull,X0,options,testContrast,datavals,[],'least_square',[1 1 1 0],[10 100 1 0],[1 1 0 0],[-20 0 0 0]);
    %         X=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 1],[2 100 1 0.2],[1 1 0 0],[-20 0 0 0]);
    X=fminsearch(@fit_weibull,X1,options,testContrast,datavals,[],'mle',[1 1 1 0],[10 100 1 1],[1 1 0 0],[-10 0 0 0]);
    fitted_yvals=1-X(4)-X(3).*exp(-(testContrast./X(2)).^X(1));
    maxRate(1,i)=1-X(4);
    minRate(1,i)=1-X(4)-X(3);
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
%         slopeNeuro(1,i)=(1-X(4)-X(3).*exp(-(30.001./X(2)).^X(1))-1-X(4)-X(3).*exp(-(29.999./X(2)).^X(1)))/0.002;
    end
    
    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
    xvalsFine=testContrast(1):0.01:testContrast(end);
    yvalsFine=1-X(4)-X(3).*exp(-(xvalsFine./X(2)).^X(1));
    if useISI==0
        if max(yvals)<0.5%out of range
            c50(1,i)=100;
        elseif min(yvals)>0.5
            c50(1,i)=0;
        else
            diffTemp=yvalsFine-0.5;
            [tempVal columnInd]=min(abs(diffTemp));
            c50(1,i)=xvalsFine(columnInd);
            %                 c50(1,i)=real(X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1)));
        end
        diffTemp=[];
        if slopeNeuro(1,i)>=0&&max(yvals)>=0.82%stimulus-evoked excitation
            diffTemp=yvalsFine-0.82;%single neurometric threshold
        elseif slopeNeuro(1,i)<0&&min(yvals)<=1-0.82%stimulus-evoked suppression
            diffTemp=yvalsFine-(1-0.82);%single neurometric threshold
        end
        if ~isempty(diffTemp)
            if max(yvals)>=0.82&&min(yvals)<=0.82%||max(yvals)>=1-0.82&&min(yvals)<=1-0.82
                [tempVal columnInd]=min(abs(diffTemp));
                threshold82higher(1,i)=xvalsFine(columnInd);
            else
                threshold82higher(1,i)=NaN;
            end
        else
            threshold82higher(1,i)=NaN;
        end
        diffTemp=[];
        if slopeNeuro(1,i)>=0&&min(yvals)<=0.18%stimulus-evoked excitation
            diffTemp=yvalsFine-0.18;%single neurometric threshold
        elseif slopeNeuro(1,i)<0&&max(yvals)>=0.82%stimulus-evoked suppression
            diffTemp=yvalsFine-0.82;%single neurometric threshold
        end
        if ~isempty(diffTemp)
            if max(yvals)>=0.18&&min(yvals)<=0.18
                [tempVal columnInd]=min(abs(diffTemp));
                threshold82lower(1,i)=xvalsFine(columnInd);
            else
                threshold82lower(1,i)=NaN;
            end
        else
            threshold82lower(1,i)=NaN;
        end
    elseif useISI==1%read C50/PNE as being point where AUROC value is 0.75, instead of 0.5
        if max(yvals)<0.75%out of range
            c50(1,i)=100;
        elseif min(yvals)>0.75
            c50(1,i)=0;
        else
            diffTemp=yvalsFine-0.75;
            [tempVal columnInd]=min(abs(diffTemp));
            c50(1,i)=xvalsFine(columnInd);
            %                 c50(1,i)=real(X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1)));
        end
        diffTemp=[];
        if slopeNeuro(1,i)>=0&&max(yvals)>=0.82%stimulus-evoked excitation
            diffTemp=yvalsFine-0.82;%single neurometric threshold
        elseif slopeNeuro(1,i)<0&&min(yvals)<=1-0.82%stimulus-evoked suppression
            diffTemp=yvalsFine-(1-0.82);%single neurometric threshold
        end
        if ~isempty(diffTemp)
            if max(yvals)>=0.82&&min(yvals)<=0.82%||max(yvals)>=1-0.82&&min(yvals)<=1-0.82
                [tempVal columnInd]=min(abs(diffTemp));
                threshold82higher(1,i)=xvalsFine(columnInd);
            else
                threshold82higher(1,i)=NaN;
            end
        else
            threshold82higher(1,i)=NaN;
        end
        diffTemp=[];
        if slopeNeuro(1,i)>=0&&min(yvals)<=0.18%stimulus-evoked excitation
            diffTemp=yvalsFine-0.18;%single neurometric threshold
        elseif slopeNeuro(1,i)<0&&max(yvals)>=0.82%stimulus-evoked suppression
            diffTemp=yvalsFine-0.82;%single neurometric threshold
        end
        if ~isempty(diffTemp)
            if max(yvals)>=0.18&&min(yvals)<=0.18
                [tempVal columnInd]=min(abs(diffTemp));
                threshold82lower(1,i)=xvalsFine(columnInd);
            else
                threshold82lower(1,i)=NaN;
            end
        else
            threshold82lower(1,i)=NaN;
        end
    end
end
hold on
if plotDiffC50_30==1
    diffc50(1,i)=abs(c50(1,i)-30);
else
    diffc50=[];
end
if strcmp(ROCanalysisType,'old')
%     line(c50(1,i),0:0.01:1,'Color','b');
    plot([c50(1,i) c50(1,i)],[0 1],'Color','b');
else
    if useISI==0
        line(c50(1,i),0:0.01:1,'Color','r');
        plot([c50(1,i) c50(1,i)],[0 1],'Color','r');
    elseif useISI==1
%         line(threshold82higher(1,i),0:0.01:1,'Color','r');
        plot([threshold82higher(1,i) threshold82higher(1,i)],[0 1],'Color','r');
    end
end
if alex_fit==1
    yvals=1-exp(-((xvals/X(1)).^X(2)));
else
    yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
end
%         yvals=(yvals-.005)/.99;
%         yvals=log(yvals);
if strcmp(ROCanalysisType,'old')
    plot(xvals,yvals,'Color','b');
elseif strcmp(ROCanalysisType,'new')
    plot(xvals,yvals,'Color','r');
else
    plot(xvals,yvals,'Color','k');
end
xlim([0 testContrast(end)+10]);