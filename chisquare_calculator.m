function chisq=chisquare_calculator(y1,y2)
%Written by Xing 17/02/14 to calculate maximum likelihood estimate value.
%y1 contains observed values, y2 contains fitted values.
chisq=0;
for j=1:length(y1)
    chisq=chisq+((y1(j)-y2(j))^2)/y2(j);%calculate Chi-Square goodness of fit statistic
    %as sum of [(O-E)^2/variance], where variance is simply set to 1,
    %under the assumption that it is constant across data points.
end

