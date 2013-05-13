 function err=n_r(X,x,y)

 Rmax=X(1);
 c50=X(2);
 n=X(3);%slope
 offset=X(4);
 
 err=sum((y-(Rmax*(x.^n./(x.^n+c50^n))+offset)).^2);
 
% if Rmax>lim(1)||c50>lim(2)||c50<0.1||n<0.01||n>90||offset<lim(2)
if Rmax>max(y)||c50>max(x)||c50<0.1||n<-30||n>30%||offset<min(y)
     err=100000;
 end

 