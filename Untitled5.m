
W_glider_lisse = smoothdata(explorer_all.W_glider_filter,'SmoothingFactor',0.02);
Pitch_lisse = smoothdata(explorer_all.pitch_filter','SmoothingFactor',0.03);
Temp_lisse = smoothdata(explorer_all.temp,'SmoothingFactor',0.02);


figure()
plot(explorer_all.W_glider_filter)
hold on 
plot(W_glider_lisse,'-+')
legend('data','lisse')

figure()
plot(explorer_all.pitch_filter')
hold on 
plot(Pitch_lisse,'-+')
legend('data','lisse')

figure()
plot(explorer_all.temp)
hold on 
plot(Temp_lisse,'-+')
legend('data','lisse')

