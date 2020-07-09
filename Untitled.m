figure()
subplot(3,1,1)
plot(Triplet(1,:),'Color','#D95319','LineWidth',2)
title('V0 ')
ylim([0.05937 0.05951])
subplot(3,1,2)
plot(Triplet(2,:),'Color','#77AC30','LineWidth',2)
title('eps')
ylim([3.5e-10 5e-10])
subplot(3,1,3)
plot(Triplet(3,:),'Color','#0072BD','LineWidth',2)
title('Cd')
ylim([0.08 0.17])