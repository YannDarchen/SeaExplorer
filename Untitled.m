figure('Name','desc\_1')
subplot(4,1,1)
plot(Quintuplet_desc_1(1,:),'Color','#D95319','LineWidth',2)
title('V0')
subplot(4,1,2)
 plot(Quintuplet_desc_1(2,:),'Color','#77AC30','LineWidth',2)
 title('eps')

subplot(4,1,3)
plot(Quintuplet_desc_1(3,:),'Color','#0072BD','LineWidth',2)
title('Cd')
  subplot(4,1,4)
   plot(Fb_desc_1,'Color','#77AC30','LineWidth',2)
   title('Precision')

figure('Name','desc\_2')
subplot(3,1,1)
plot(Quintuplet_desc_2(1,:),'Color','#D95319','LineWidth',2)
title('V0')
subplot(3,1,2)
plot(Quintuplet_desc_2(2,:),'Color','#77AC30','LineWidth',2)
title('eps')
subplot(3,1,3)
plot(Quintuplet_desc_2(3,:),'Color','#0072BD','LineWidth',2)
title('Cd')

figure('Name','mont\_1')
subplot(3,1,1)
plot(Quintuplet_mont_1(1,:),'Color','#D95319','LineWidth',2)
title('V0')
subplot(3,1,2)
plot(Quintuplet_mont_1(2,:),'Color','#77AC30','LineWidth',2)
title('eps')
subplot(3,1,3)
plot(Quintuplet_mont_1(3,:),'Color','#0072BD','LineWidth',2)
title('Cd')

figure('Name','mont\_2')
subplot(3,1,1)
plot(Quintuplet_mont_2(1,:),'Color','#D95319','LineWidth',2)
title('V0')
subplot(3,1,2)
plot(Quintuplet_mont_2(2,:),'Color','#77AC30','LineWidth',2)
title('eps')
subplot(3,1,3)
plot(Quintuplet_mont_2(3,:),'Color','#0072BD','LineWidth',2)
title('Cd')
