tab_pitch1=tableau_pitch(:,1);
to_ign=isnan(tab_pitch1);
tab_pitch1(to_ign ~= 0)=[];
tab_pitch1_1=detrend(tab_pitch1)-25;
figure
plot(tab_pitch1,-pi)
hold on 
plot(tab_pitch1_1,-pi)