 function err=weibull_pointfive_one(p,x,response,var_max_limit,max_limit,var_min_limit,min_limit,err_type)
 %Fit curve to neurometric or psychometric data for comparison of 63%
 %thresholds
 global ERR
 %err=	sqrt	(sum((y-(1-.5*exp(-((x/vars(1)).^vars(2))))).^2));		  
 %err=	sum(sqrt( (y-(1-.5*exp(-( (x/vars(1)).^vars(2)) ) ) ) .^2));
 err=	sum(((response-(1-.5*exp(-((x/p(1)).^p(2)))))).^2);	
 
% -------------------------------------------------------------------------
 ERR=err;
 if strcmp(err_type,'mle')==1
    y=1-.5.*exp(-(conditions./p(2)).^p(1));
    y=y*.99+.005;
    err=sum(response.*log(y) + (1-response).*log(1-y));
 end

 % limit % -----------------------------------------------------------------
if var_max_limit(1) ~= 0
    if p(1) > max_limit(1)
        err = 10^1000;
    end
end
if var_max_limit(2) ~= 0
    if p(2) > max_limit(2)
        err = 10^1000;
    end
end
% ---------------------------
if var_min_limit(1) ~= 0
    if p(1) < min_limit(1)
        err = 10^1000;
    end
end
if var_min_limit(2) ~= 0
    if p(2) < min_limit(2)
        err = 10^1000;
    end
end