 function err=weibull_pointfive_one(p,x,y,var_max_limit,max_limit,var_min_limit,min_limit)
 %Fit curve to neurometric or psychometric data for comparison of 63%
 %thresholds
 global ERR
 %err=	sqrt	(sum((y-(1-.5*exp(-((x/vars(1)).^vars(2))))).^2));		  
 %err=	sum(sqrt( (y-(1-.5*exp(-( (x/vars(1)).^vars(2)) ) ) ) .^2));
 err=	sum(((y-(1-.5*exp(-((x/p(1)).^p(2)))))).^2);	
 ERR=err;
 
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