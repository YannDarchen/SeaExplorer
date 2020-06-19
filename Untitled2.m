figure()
L=ones(length(tab_descent_1.tab(:,8)))*(-0.005);
U=ones(length(tab_descent_1.tab(:,8)))*(0.005);
plot(tab_descent_1.tab(:,8),'-+')
hold on 
plot(L,'--r')
plot(U,'--r')