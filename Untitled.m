t = 0:20;
x = 3*sin(t) + t;
y = detrend(x);
plot(t,x,t,y,t,x-y,':k')
legend('Input Data','Detrended Data','Trend','Location','northwest') 