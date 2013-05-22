function [slopeNeuro,c50,diffc50,minRate,maxRate,chSSE]=weibull_fitting(datavals,sampleContrast,testContrast,ROCanalysisType,i,slopeNeuro,chSSE,c50,minRate,maxRate,diffc50,plotDiffC50_30,calculateTangent)
xvals=testContrast(1):1:testContrast(end);
if sum(datavals(1:3))<=sum(datavals(end-2:end))
    X0=[2 30 0.2 0.1];
elseif sum(datavals(1:3))>sum(datavals(end-2:end))
    X0=[-2 30 0.2 0.1];
end
options = optimset('Display','off','MaxFunEvals',10^4,'MaxIter',10^4,'TolFun',1.0E-6,'TolX',1.0E-6);
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
end
c50(1,i)=real(X(2).*(-log((0.5-X(4))/X(3))).^(1/X(1)));
hold on
if plotDiffC50_30==1
    diffc50(1,i)=abs(c50(1,i)-30);
else 
    diffc50=[];
end
if strcmp(ROCanalysisType,'old')
    line(c50(1,i),0:0.01:1,'Color','b');
else
    line(c50(1,i),0:0.01:1,'Color','r');
end
yvals=1-X(4)-X(3).*exp(-(xvals./X(2)).^X(1));
%         yvals=(yvals-.005)/.99;
%         yvals=log(yvals);
if strcmp(ROCanalysisType,'old')
    plot(xvals,yvals,'Color','b');
elseif strcmp(ROCanalysisType,'new')
    plot(xvals,yvals,'Color','r');
else
    plot(xvals,yvals,'Color','k');
end