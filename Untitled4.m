figure('Name','Wglider and Wmodel all','NumberTitle','off','Units','normalized','Position',[0,0.05,0.33,0.38])
 % plot(explorer_all.time, explorer_all.W_glider_filter,'LineWidth',1.5)
 %yyaxis left 
 plot(explorer_all.time, explorer_all.W_glider_filter,'LineWidth',1.5)
  title('Vitesses verticales du glider all')
  hold on 
  plot(explorer_all.time,W_model_all,'LineWidth',1.5)
  %yyaxis right
  plot(explorer_all.time,explorer_all.pitch_filter.*(pi/180),'-+','LineWidth',1)
  xticks(explorer_all.time(1):0.2:explorer_all.time(end))
  xticklabels(datestr(explorer_all.time(1):0.2:explorer_all.time(end)))
  %datetick('x',0,'keepticks')
  legend('W\_glider','W\_model','pitch')