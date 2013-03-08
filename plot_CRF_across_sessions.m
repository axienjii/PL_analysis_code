function [slopeNeuro,c50,diffc50,minRate,maxRate,chSSE,yLimCRF]=plot_CRF_across_sessions(i,numsessions,testContrast,crfvals,calculateTangent,plotDiffC50_30,slopeNeuro,c50,diffc50,minRate,maxRate,chSSE,yLimCRF)

xvals=testContrast(1):1:testContrast(end);
subplot(ceil(numsessions/5),5,i);
fittingFunction='n_r';
options = optimset('Display','off','MaxFunEvals',10^6,'MaxIter',10^6,'TolFun',1.0E-6,'TolX',1.0E-6);
plot(testContrast,crfvals,'ok');
hold on
X0=[max(crfvals) 30 0.5 min(crfvals)];
if mean(crfvals(1:3))>mean(crfvals(12:14))%||chNum==13.2||chNum==24||chNum==42
    X0=[max(crfvals) 30 -0.1 min(crfvals)];%negative slope
end
[X]=fminsearch(fittingFunction,X0,options,testContrast,crfvals);

if calculateTangent==0
    slopeNeuro(1,i)=X(3);
elseif calculateTangent==1
    slopeNeuro(1,i)=X(2)^X(3)*30^(X(3)-1)/(30^X(3)+X(2)^X(3))^2;
end
c50(1,i)=X(2);
minRate(1,i)=X(4);
maxRate(1,i)=X(1);
%     [X,fval]=fminsearch('weib_sim_min_max',X0,options,testContrast,crfvals);
%     slopeNeuro(1,i)=X(2);
%     c50(1,i)=X(1).*(-log(0.5)).^(1/X(2));
%     yvals=max(crfvals)-(max(crfvals)-min(crfvals)).*exp(-((xvals./X(1
%     )).^X(2)));
if plotDiffC50_30==1
    diffc50(1,i)=abs(c50(1,i)-30);
else
    diffc50=[];
end
yvals=X(1)*(xvals.^X(3)./(xvals.^X(3)+X(2).^X(3)))+X(4);
fitted_yvals=X(1)*(testContrast.^X(3)./(testContrast.^X(3)+X(2).^X(3)))+X(4);
residuals = crfvals-fitted_yvals;
sseCRF=sum(residuals.^2);
chSSE(i,:)=[1:numsessions sseCRF];
plot(xvals,yvals,'r');
%         ylim([0,max(crfvals)]);
% ylim([0 1]);
xlim([0 max(testContrast)]);
subplottitle=num2str(i);
title(subplottitle);
[yLimCRF(i,1:2)]=get(gca,'YLim');