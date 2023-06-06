[T,Res]=lyapunov(3,@lorenz_ext,@ode45,0,0.1,40000,[0.1 0.1 0.1],10);
plot(T,Res);
title('Dynamics of Lyapunov exponents');
xlabel('Time'); ylabel('Lyapunov exponents');

