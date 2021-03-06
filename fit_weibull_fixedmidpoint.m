function [err] = fit_weibull_fixedmidpoint(p,conditions,response,standard_dev,err_type,var_max_limit,max_limit,var_min_limit,min_limit,minMax)

% conditions = [0:10:90];
% response
% standard_dev or []
% error_type: 'least_square' 
% var_max_limit= [1 1 0]
% max_limit= [6 1 0]
% var_min_limit= [1 1 0]
% min_limit= [5 0.8 0]

offset = response - (minMax(2)-minMax(1).*exp(-(conditions./30).^p(1)));

% p(1) : slope
% p(2) : max
% p(3) : min

err = 0;
% -------------------------------------------------------------------------
if (strcmp(err_type,'least_square') == 1)    
    if isempty(standard_dev) == 0
        for i = 1 : length(conditions)
            if standard_dev(i) ~= 0
                err = err + (offset(i)^2);
            end
        end
    end    
    if isempty(standard_dev) == 1
        err = sum(offset.^2);
    end    
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
if strcmp(err_type,'mle')==1
    y=(1-p(4))-p(3).*exp(-(conditions./p(2)).^p(1));
%     y=y*.99+.005;
%     err=sum(response.*log(y) + (1-response).*log(1-y));
    err=-prod( y.^response .* (1-y).^(1-response));
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
if (strcmp(err_type,'chi2') == 1)
    if isempty(standard_dev) == 0
        for i = 1 : length(conditions)
            if (standard_dev(i) < 1)
                standard_dev(i) = 1; 
            end        
        err = err + ((offset(i)^2)/(standard_dev(i)^2));
        end
    end
    if isempty(standard_dev) == 1
        error('you need standard deviation for minimizing chi square error')
    end
end
% -------------------------------------------------------------------------

% limit % -----------------------------------------------------------------
if var_max_limit(1) ~= 0
    if p(1) > max_limit(1)
        err = 10^1000;
    end
end
% ---------------------------
if var_min_limit(1) ~= 0
    if p(1) < min_limit(1)
        err = 10^1000;
    end
end