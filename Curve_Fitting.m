%% Curve fitting for a non-linear equation
clear
clc
clf
data = table2array(readtable('data-assignment04.csv','Range','A1:C121'));
n = 20;  % Number of data points to be sampled
idx = randperm(121);  % Non-repeating random numbers from 1 to 121 are generated
% A heuristic method is used to collect sample points from the domain such
% that L2 norm of error is less for any random point selection.
for i =1:n
   if idx(i) < 30
       p = idx(i)+30;
   elseif idx(i) > 90
       p = idx(i)-30;
   else 
       p = idx(i);  
   end
   q(i) = p;
    w = unique(q);
    z = abs(numel(q)-numel(w));
    if  z > 0
        r = [w randperm(121,z)];
        q = r;       
end
end

for j = 1:n  % Allocation of random data points 
    xyz(j,:) = data(q(j),:); 
end
   
x1 = xyz(:,1); x2 = xyz(:,2); y = xyz(:,3);
x = [x1 x2];
% Function Definition
fun = @(a,x)a(1)*x(:,1).^2 + a(2)*x(:,2).^2 + a(3)*x(:,1) + a(4)*x(:,2) + a(5);
b0 = [1 1 1 1 1];       % Coefficients Beta_i are initialized
[beta,resnorm,L2,exitflag,output] = lsqcurvefit(fun,b0,x,y);  %trust region algortithm is used to
                                                             %value of Beta
y_ = beta(1)*x(:,1).^2 + beta(2)*x(:,2).^2 + beta(3)*x(:,1) + beta(4)*x(:,2) + beta(5);
L2 = norm(L2,2)     % L2 norm of the residue is calculated
figure
scatter3(x1,x2,y);
hold on             % 3D plot from observed values
scatter3(x1,x2,y_); % 3D plot from calculated values
title("Curve fitting for "+numel(unique(q))+" points")
hold off
legend('Given Data','Curve-fitted data')
xlabel('x1')
ylabel('x2')
zlabel('y')


