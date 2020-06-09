% TF = isoutlier(explorer.W_glider,'median');
% explorer.W_glider_filter = explorer.W_glider; 
% explorer.W_glider_filter(TF)=NaN;
% explorer.W_glider_filter=fillmissing(explorer.W_glider_filter,'next');
% figure()
% hold on 
% plot(explorer.time(1:end-5),explorer.W_glider,'-+')
% plot(explorer.time(1:end-5),explorer.W_glider_filter,'-+')
% legend('data','filter')
ind_neg = explorer.W_glider < 0;
negatif = explorer.W_glider(ind_neg);


% t = datetime(2017,1,1,0,0,0) + hours(0:length(negatif)-1);
% TF = isoutlier(negatif,'movmedian',hours(500),'SamplePoints',t);
% figure()
% hold on 
% plot(t,negatif,'-+')
% plot(t(TF),negatif(TF),'x')
% legend('Data','Outlier')

b=length(negatif);
[TF,L,U,C] = isoutlier(negatif,'percentiles',[3 100]);

figure()
plot(explorer.time(1:b),negatif,'-+',explorer.time(TF),explorer.W_glider(TF),'x',explorer.time(1:b),L*ones(1,b),explorer.time(1:b),U*ones(1,b),explorer.time(1:b),C*ones(1,b))
legend('Original Data','Outlier','Lower Threshold','Upper Threshold','Center Value')
explorer.W_glider(TF)=NaN;
explorer.W_glider = fillmissing(explorer.W_glider,'next');
figure()
plot(explorer.time(1:end-5),explorer.W_glider)


