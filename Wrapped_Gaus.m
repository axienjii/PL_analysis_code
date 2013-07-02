function [Err] = Wrapped_Gaus(p,Response,StandDev)
%input arg 'p':
%p(1) difference between max and min response
%p(2) orientation at which max response obtained
%p(3) width at HWHH
%p(4) min response

x = [0:15:165];

Offset = Response-(p(1).*(exp((-(x-p(2)+180*(-5)).^2)./(2*p(3).^2))+exp((-(x-p(2)+180*(-4)).^2)./(2*p(3).^2))+...
         exp((-(x-p(2)+180*(-3)).^2)./(2*p(3).^2))+exp((-(x-p(2)+180*(-2)).^2)./(2*p(3).^2))+...
         exp((-(x-p(2)+180*(-1)).^2)./(2*p(3).^2))+exp((-(x-p(2)+180*(0)).^2)./(2*p(3).^2))+...
         exp((-(x-p(2)+180*(1)).^2)./(2*p(3).^2))+exp((-(x-p(2)+180*(2)).^2)./(2*p(3).^2))+...
         exp((-(x-p(2)+180*(3)).^2)./(2*p(3).^2))+exp((-(x-p(2)+180*(4)).^2)./(2*p(3).^2))+...
         exp((-(x-p(2)+180*(5)).^2)./(2*p(3).^2)))+p(4));

% chi2
% Chi2 = 0;
% for i = 1 : 12
%     if (StandDev(i) < 1)
%        StandDev(i) = 1; 
%     end        
%     Chi2 = Chi2 + ((Offset(i)^2)/(StandDev(i)^2));
% 
% end
% Err = Chi2;

% if p(1)>1.5*(max(Response)-min(Response))
%     Err=10^10;
% else    
%     Err = Chi2;
% end

% % LS_1
% Err = sum(Offset.^2);

% LS_2
Err = 0;
for i = 1 : 12
    if StandDev(i) ~= 0
        Err = Err + (Offset(i)^2);
    end
end

% ABS
% Err = sum(abs(Offset));