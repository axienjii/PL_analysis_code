function check_min_max_weibull
xvals=0:1:100;
figure
yvals=1-0-1.*exp(-(xvals./30).^X(1));
plot(xvals,yvals,'k-','Marker','none');
hold on
yvals=1-0-0.5.*exp(-(xvals./30).^X(1));
plot(xvals,yvals,'r-','Marker','none');
hold on
yvals=1-0-0.3.*exp(-(xvals./30).^X(1));
plot(xvals,yvals,'k-','Marker','none');
hold on
yvals=1-0.5-0.3.*exp(-(xvals./30).^X(1));
plot(xvals,yvals,'g-','Marker','none');
hold on
yvals=1-0.6-0.3.*exp(-(xvals./30).^X(1));
plot(xvals,yvals,'r-','Marker','none');
hold on
yvals=1-0.5-0.2.*exp(-(xvals./30).^X(1));
plot(xvals,yvals,'b-','Marker','none');
hold on
yvals=1-0.4-0.2.*exp(-(xvals./30).^X(1));
plot(xvals,yvals,'k-','Marker','none');
hold on
yvals=1-0.0-0.2.*exp(-(xvals./30).^X(1));
plot(xvals,yvals,'g-','Marker','none');