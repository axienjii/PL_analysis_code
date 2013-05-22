function exampleFig_gauss(cond1t,cond1s,cond14t,cond14s)
maxRate=max([cond1t cond1s cond14t cond14s]);
cond1st=[cond1s cond14s];
numBins=20;
bins=0:ceil(maxRate)/numBins:ceil(maxRate);
for binInd=1:length(bins)-1
    higher=cond1t>bins(binInd);
    lower=cond1t<=bins(binInd+1);
    cond1tBins(binInd)=sum(higher+lower==2);
    higher=cond14t>bins(binInd);
    lower=cond14t<=bins(binInd+1);
    cond14tBins(binInd)=sum(higher+lower==2);
    higher=cond1st>bins(binInd);
    lower=cond1st<=bins(binInd+1);
    condstBins(binInd)=sum(higher+lower==2);
end
sigma=3;
mu=0;
vector=cond1tBins;
fit_val=[(-3*sigma):1:(3*sigma)];  
y=ones(1, length(fit_val));             
y=(y*(1/(sigma*sqrt(2*pi))).*exp(-(((fit_val-mu).^2)/(2*sigma*sigma))));
fitted_vector=filter2(y,vector);
midBins=bins(2:end)-ceil(maxRate)/(numBins*2);
plot(midBins,fitted_vector,'ro','MarkerSize',5);hold on
[a b]=max(fitted_vector);
plot([midBins(b) midBins(b)],[0 a],'r','LineWidth',2);
vector=cond14tBins;
fit_val=[(-3*sigma):1:(3*sigma)];  
y=ones(1, length(fit_val));             
y=(y*(1/(sigma*sqrt(2*pi))).*exp(-(((fit_val-mu).^2)/(2*sigma*sigma))));
fitted_vector=filter2(y,vector);
midBins=bins(2:end)-ceil(maxRate)/(numBins*2);
plot(midBins,fitted_vector,'bo','MarkerSize',5);hold on
[a b]=max(fitted_vector);
plot([midBins(b) midBins(b)],[0 a],'b','LineWidth',2);
vector=condstBins;
fit_val=[(-3*sigma):1:(3*sigma)];  
y=ones(1, length(fit_val));             
y=(y*(1/(sigma*sqrt(2*pi))).*exp(-(((fit_val-mu).^2)/(2*sigma*sigma))));
fitted_vector=filter2(y,vector);
midBins=bins(2:end)-ceil(maxRate)/(numBins*2);
plot(midBins,fitted_vector,'ko','MarkerSize',5);hold on
[a b]=max(fitted_vector);
plot([midBins(b) midBins(b)],[0 a],'k','LineWidth',2);