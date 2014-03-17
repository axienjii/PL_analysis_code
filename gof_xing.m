function chiSq=gof_xing(x,xref)
%Written by Xing 10/02/14.
%Calculates goodness of fit between fitted values and actual data, based on
%mean square error.
chiSq=0;
for count=1:length(x)
    chiSq=chiSq+((x(count)-xref(count))^2)/xref(count);
end
% chiSqStat=chiSq/length(x);