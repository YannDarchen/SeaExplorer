close all;
figure('Name','desc\_1')

subplot(3,1,1)
plot(Quintuplet_desc_1(1,:),'Color','#D95319','LineWidth',2)
title('V0  descent\_1')
%ylim([0.0588 0.0591])
subplot(3,1,2)
plot(Quintuplet_desc_1(2,:),'Color','#77AC30','LineWidth',2)
title('eps')
ylim([3.5e-10 5.465e-10])
subplot(3,1,3)
plot(Quintuplet_desc_1(3,:),'Color','#0072BD','LineWidth',2)
title('Cd')
ylim([0.04 0.15])
%  subplot(4,1,4)
%  plot(Fb_desc_1,'Color','#77AC30','LineWidth',2)
%  title('Precision')

figure('Name','desc\_2')
subplot(3,1,1)
plot(Quintuplet_desc_2(1,:),'Color','#D95319','LineWidth',2)
title('V0  descent\_2')
%ylim([0.0588 0.0591])
subplot(3,1,2)
plot(Quintuplet_desc_2(2,:),'Color','#77AC30','LineWidth',2)
title('eps')
ylim([3.5e-10 5.465e-10])
subplot(3,1,3)
plot(Quintuplet_desc_2(3,:),'Color','#0072BD','LineWidth',2)
title('Cd')
ylim([0.04 0.15])

figure('Name','mont\_1')
subplot(3,1,1)
plot(Quintuplet_mont_1(1,:),'Color','#D95319','LineWidth',2)
title('V0   montee\_1')
%ylim([0.0588 0.0591])
subplot(3,1,2)
plot(Quintuplet_mont_1(2,:),'Color','#77AC30','LineWidth',2)
title('eps')
ylim([3.5e-10 5.465e-10])
subplot(3,1,3)
plot(Quintuplet_mont_1(3,:),'Color','#0072BD','LineWidth',2)
title('Cd')
ylim([0.04 0.15])

figure('Name','mont\_2')
subplot(3,1,1)
plot(Quintuplet_mont_2(1,:),'Color','#D95319','LineWidth',2)
title('V0   montee\_2')
%ylim([0.0588 0.0591])
subplot(3,1,2)
plot(Quintuplet_mont_2(2,:),'Color','#77AC30','LineWidth',2)
title('eps')
ylim([3.5e-10 5.465e-10])
subplot(3,1,3)
plot(Quintuplet_mont_2(3,:),'Color','#0072BD','LineWidth',2)
title('Cd')
ylim([0.04 0.15])

figure('Name','desc\_1')

subplot(3,1,1)
plot(Quintuplet_desc_1(1,:),'Color','#D95319','LineWidth',2)
title('V0  descent\_1')
%ylim([0.0588 0.0591])
subplot(3,1,2)
plot(Quintuplet_desc_1(4,:),'Color','#77AC30','LineWidth',2)
title('aw')
ylim([2 3])
subplot(3,1,3)
plot(Quintuplet_desc_1(5,:),'Color','#0072BD','LineWidth',2)
title('cdw')
ylim([2 3])

figure('Name','mont\_1')

subplot(3,1,1)
plot(Quintuplet_mont_1(1,:),'Color','#D95319','LineWidth',2)
title('V0  montee\_1')
%ylim([0.0588 0.0591])
subplot(3,1,2)
plot(Quintuplet_mont_1(4,:),'Color','#77AC30','LineWidth',2)
title('aw')
ylim([2 3])
subplot(3,1,3)
plot(Quintuplet_mont_1(5,:),'Color','#0072BD','LineWidth',2)
title('cdw')
ylim([2 3])
