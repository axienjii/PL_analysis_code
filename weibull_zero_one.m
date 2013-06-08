 function err=weibull_zero_one(vars,x,y)
 
 global ERR
 %err=	sqrt	(sum((y-(1-.5*exp(-((x/vars(1)).^vars(2))))).^2));		  
 %err=	sum(sqrt( (y-(1-.5*exp(-( (x/vars(1)).^vars(2)) ) ) ) .^2));		  
 err=	sum(((y-(1-exp(-((x/vars(1)).^vars(2)))))).^2);	
 ERR=err;