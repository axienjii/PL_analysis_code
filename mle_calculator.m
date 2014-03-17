function xmle=mle_calculator(n,y1,y2)
%Written by Xing 17/02/14 to calculate maximum likelihood estimate value.
%y1 contains observed values, y2 contains fitted values.
xmle=0;
for j=1:length(y1)
    xmle=xmle+n(j)*(y1(j)*log(y2(j))+(1-y1(j))*log(1-y2(j)));%use natural log
        %equivalent to xmle=xmle+n(j)*(log(obsp^obsp)+log(obsq^obsq)),
        %which is equivalent to
        %xmle=xmle+n(j)*((log(obsp))^obsp)+(log(obsq))^obsq))
        %Note that it doesn't matter whether the 'to the power' is applied to the
        %whole log term or to just the quantity within the brackets of the
        %log term, i.e. log(2^3) == (log(2))^3 == 3*log(2)
        %Also note that 'log' function in matlab is natural log, whereas 'log10' is log with base 10
end

