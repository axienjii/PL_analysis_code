 function err=weibull(vars,x,y)
 
 global ERR
 %err=	sqrt	(sum((y-(1-.5*exp(-((x/vars(1)).^vars(2))))).^2));		  
 %err=	sum(sqrt( (y-(1-.5*exp(-( (x/vars(1)).^vars(2)) ) ) ) .^2));		  
 err=	sum(((y-(1-.5*exp(-((x/vars(1)).^vars(2)))))).^2);	
 if vars(1)<0.1 | vars(2)<0.01
     err=10^100;
 end
     
 ERR=err;