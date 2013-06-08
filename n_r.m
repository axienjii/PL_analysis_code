 function err=n_r(X,x,y,var_max_limit,max_limit,var_min_limit,min_limit)

 Rmax=X(1);
 c50=X(2);
 n=X(3);%slope
 offset=X(4);
 
 err=sum((y-(Rmax*(x.^n./(x.^n+c50^n))+offset)).^2);
 
% if Rmax>lim(1)||c50>lim(2)||c50<0.1||n<0.01||n>90||offset<lim(2)
if Rmax>max(y)||c50>max(x)||c50<0.1||n<-30||n>30%||offset<min(y)
     err=100000;
end

if nargin>3
    % limit % -----------------------------------------------------------------
    if var_max_limit(1) ~= 0
        if X(1) > max_limit(1)
            err = 10^1000;
        end
    end
    if var_max_limit(2) ~= 0
        if X(2) > max_limit(2)
            err = 10^1000;
        end
    end
    if var_max_limit(3) ~= 0
        if X(3) > max_limit(3)
            err = 10^1000;
        end
    end
    if var_max_limit(4) ~= 0
        if X(4) > max_limit(4)
            err = 10^1000;
        end
    end
    % ---------------------------
    if var_min_limit(1) ~= 0
        if X(1) < min_limit(1)
            err = 10^1000;
        end
    end
    if var_min_limit(2) ~= 0
        if X(2) < min_limit(2)
            err = 10^1000;
        end
    end
    if var_min_limit(3) ~= 0
        if X(3) < min_limit(3)
            err = 10^1000;
        end
    end
    if var_min_limit(4) ~= 0
        if X(4) < min_limit(4)
            err = 10^1000;
        end
    end
end