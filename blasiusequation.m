% Blasius equation over Flat Plate
clc 
clear
Vinf = 5;
visc = 18.6e-6;
f = zeros(4,1);
k = 0;
while k >= 0
[eta,F] = blas(f,k);
if round(F(end,2),2) == 1
    figure
    plot(F(:,2),eta);
    xlabel("f'");
    ylabel('eta');
    break
else
    k = k+0.01;
end
end
x = linspace(0,1,77);
X = meshgrid(x);

y = eta*sqrt(x*visc/Vinf);

tsi = F(:,1)*sqrt(x*Vinf*visc);
figure
plot(x,tsi(1,:),'b');
hold on 
plot(x,tsi(30,:),'r');
hold on 
plot(x,tsi(50,:),'k');
hold on 
plot(x,tsi(77,:),'c');
hold on 

del = 4.9*sqrt(visc/Vinf)*sqrt(x);
plot(x,del,'y');
legend('Streamline 1','Streamline 2','Streamline 3','Streamline 4','Boundary layer');
title('Streamline variation');
hold off
[Vy,Vx] = gradient(tsi);
Vy = -Vy;
figure
quiver(X,y,Vx,Vy);
xlabel('x axis (along plate)');
ylabel('y axis');
title('Velocity profile along flat plate');


function [eta,F] = blas(~,y)

    df = @(e,f)[f(2);f(3); (-f(1).*f(3))./2];
    [eta,F] = ode45(df, [0 8], [0 0 y]);
end   