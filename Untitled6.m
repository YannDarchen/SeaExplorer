figure()
plot(explorer_all.M.*explorer_all.W_glider_acc,'-+')
hold on 
L=ones(length(explorer_all.time),1)*0.007;
K=ones(length(explorer_all.time),1)*(-0.007);
plot(L,'--r')
plot(K,'--r')
title('M*acceleration')
legend('acceleration','pitch')
