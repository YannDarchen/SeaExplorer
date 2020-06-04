close all;clc;clear;clear global;
t = 1:0.2:15;
A = sin(2*pi*t) + cos(2*pi*0.5*t);
Anoise = A + 0.5*rand(1,length(t));
Anoise(36) = 20;
Amedian = smoothdata(Anoise,'movmedian');
plot(t,Anoise,t,Amedian)
axis tight
legend('Noisy Data','Moving Median')
TF = isoutlier(Anoise);
ind = find(TF)
Aoutlier = Anoise(ind)
Afill = filloutliers(Anoise,'next');
plot(t,Anoise,t,Afill)
axis tight
legend('Noisy Data with Outlier','Noisy Data with Filled Outlier')