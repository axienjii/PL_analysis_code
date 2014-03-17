function[residual, xmle, chisq]=Max_like_estimates_weibull_xing(n,y1,y2)
%edited from Alex's function, Max_like_estimates_xing_weibull on 20/01/14.
%Annotates lines of code with explanations.

%Refer to Watson (1979), "Probability sumamtion over time." Justifies and
%derives MLE equations for a Weibull context.

%Also refer to In Jae Myung (2003). "Tutorial on
%maximumum likelihood estimation." Journal of Mathematical Psychology. 

if length(n)==1
    n=ones(1,length(y1))+n;
end
xmle=0;
chisq=0;
% y1
% y2
% n
for j=1:length(y1)
    xpectp=y2(j);
    xpectq=1-xpectp;
    obsp=y1(j);
    obsq=1-obsp;
    if ~isinf(log((obsp^obsp)*(obsq^obsq)))
        xmle=xmle+n(j)*(log((obsp^obsp)*(obsq^obsq)));%equivalent to 'M0' in Watson (1979), eqn A6
        %This is the maximum of M under H0, where nothing is assumed about
        %the fit- the 'saturated model.'
        %equivalent to xmle=xmle+n(j)*(log(obsp^obsp)+log(obsq^obsq)),
        %which is equivalent to
        %xmle=xmle+n(j)*((log(obsp))^obsp)+(log(obsq))^obsq))
        %Note that it doesn't matter whether the 'to the power' is applied to the
        %whole log term or to just the quantity within the brackets of the
        %log term, i.e. log(2^3) == (log(2))^3 == 3*log(2)
        %Also note that 'log' function in matlab is natural log, whereas 'log10' is log with base 10
    end
    if ~isinf(log((xpectp^obsp)*(xpectq^obsq)))
        chisq=chisq+n(j)*(log((xpectp^obsp)*(xpectq^obsq)));%equivalent to 'M' in Watson (1979), eqn A5
        %This is the maximum of M under H1, where the fitted values are
        %being assessed.
    end
end
residual=2*(xmle-chisq);
%provides 'deviance' statistic- the MLE analog to summed square residuals
%in the ordinary least squares method