clear all
close all
format compact

a = -1.3;
b = 2;
c = -1.4;
d = -4;
x=0:1:testContrast(end);

y = 1 - d - c .* exp(- (x./b).^a);

dydx_ana =  a * c * exp( -(x./b).^a ) .* x.^(a-1) .* (1/b)^a;
dydx_num = diff(y) ./ min(diff(x));%approximate derivative- performed as a check 
x_num = x(1:end-1) + diff(x)/2;

figure
subplot(1,2,1)
plot(x,y)
ylabel('y')
subplot(1,2,2)
plot(x,dydx_ana)
hold on
plot(x_num,dydx_num,'r:')%dydx_ana should be more or less equal to dydx_num
ylabel('dy/dx')

dydx_30 =  a * c * exp( -(30/b)^a ) .* 30^(a-1) * (1/b)^a;