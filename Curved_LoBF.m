y=[0.1572,.147,.1265,.143,.0892,.0744,.0565,.0434,.0286]';
x=[.01,.009,.008,.007,.006,.005,.004,.003,.002]';
n = length(x);
F = [ones(n, 1),x,log(x)];
a = F\y;
yfit = F*a;
plot(x, y, "or", x, yfit, 'b')
yav = sum(y)/n;
St = sum( (y-yav).^2 );
Sr = sum( (y-yfit).^2 ); % residual error
R1 = sqrt((St-Sr) / St) 
R2 = corrcoef(yfit, y)
