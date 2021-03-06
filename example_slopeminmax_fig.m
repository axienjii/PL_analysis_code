function example_slopeminmax_fig
% y(2,1:14)=[0,0,0.0298,0.138,0.347,0.329,0.360,0.538,0.670,0.750,0.915,0.984,1,1];
% y(1,1:14)=[y(2,1:7)+0.2 y(2,8:end)-0.2];
% y(3,1:14)=[y(2,1:7)+0.2 y(2,8:end)];
% y(4,1:14)=[y(2,1:7) y(2,8:end)-0.2];
% y(:,7:8)=0.5;
% testContrast=[10 15 20 25 27 28 29 31 32 33 35 40 50 60];
y(2,:)=[0 0 0 0.39 0.4 0.42 0.45 0.47 1 1 1];
y(1,:)=[0.2 0.2 0.2 0.4 0.42 0.47 0.48 0.5 0.8 0.8 0.8];
y(3,:)=[0.2 0.2 0.2 0.39 0.41 0.44 0.46 0.48 0.8 0.8 1];
y(4,:)=[0 0.2 0.2 0.5 0.5 0.5 0.5 0.5 0.8 0.8 0.8];
testContrast=[10 15 20 27 28 30 32 33 50 55 60];
xvals=testContrast(1):1:testContrast(end);
xvals=0:1:100;
X0=[2 0.2 0.1];
X0=2;
options = optimset('Display','off','MaxFunEvals',10^4,'MaxIter',10^4,'TolFun',1.0E-6,'TolX',1.0E-6);
figure
subplotRemap=[3 4 2 6];
subplotCol='krbm';
for plotInd=1:4    
    X=fminsearch(@fit_weibull_fixedmidpoint,X0,options,testContrast,y(plotInd,:),[],'least_square',[1],[10],[1],[-20],[y(plotInd,end)-y(plotInd,1) y(plotInd,end)]);
    yvals(plotInd,:)=y(plotInd,end)-(y(plotInd,end)-y(plotInd,1)).*exp(-(xvals./30).^X(1));
%     yvals(plotInd,:)=1-(1-y(plotInd,end))-y(plotInd,1).*exp(-(xvals./30).^X(1));
    Xall(plotInd,:)=X;
    subplot(3,2,subplotRemap(plotInd));
    plot(xvals,yvals(plotInd,:),'k-','Marker','none','LineWidth',2);
    xlim([10 60]);
    xlim([0 100]);
    ylim([0 1]);
    c50(plotInd,1)=real(30.*(-log((0.5-(1-y(plotInd,end)))/y(plotInd,1))).^(1/X(1)));
    line(c50(plotInd,1),0:0.01:1,'Color','r');
end
Xall
hold on
yvals=1-0.4-0.3.*exp(-(xvals./30).^X(1));
plot(xvals,yvals,'y-','Marker','none');
