function [S, MatlabExp, MyExp, E] = rskew5(k)
format long
rng('shuffle')
n = 10^k; %k= 5 or 8 (try 20 for a smile :))
r =  randi([-5,5],1,10);
S = [ 0     r(1) r(2)  r(3)   r(4);...
     -r(1)   0   r(5)  r(6)   r(7);...
     -r(2) -r(5)  0    r(8)   r(9);...
     -r(3) -r(6) -r(8) 0      r(10);...
     -r(4) -r(7) -r(9) -r(10)  0    ];
MatlabExp = expm(n*S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing my exponential
A = charpoly(S); a = A(3); b = A(5);
t1 = sqrt((a-sqrt(a^2-4*b))/2);
t2 = sqrt((a+sqrt(a^2-4*b))/2);

c1 = (t2^3*sin(n*t1)-t1^3*sin(n*t2))/(sqrt(a^2*b-4*b^2));
c2 = (t2^4*(1-cos(n*t1))-t1^4*(1-cos(n*t2)))/(b*sqrt(a^2-4*b));
c3 = (t2*sin(n*t1)-t1*sin(n*t2))/(sqrt(a^2*b-4*b^2));
c4 = (t2^2*(1-cos(n*t1))-t1^2*(1-cos(n*t2)))/((b*sqrt(a^2-4*b)));

MyExp = eye(5) + c1*S + c2*S^2 + c3*S^3 + c4*S^4;
rel_error = MyExp-MatlabExp;
E = norm(rel_error);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%508projectplot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
errors = zeros(1,21);
figure(1)
xlabel('k-values'); ylabel('$\frac{1}{100}\log||M-\hat{M}||$','Interpreter','latex');
title('$k$ vs. $\frac{1}{100}\log||M-\hat{M}||$ Plot','Interpreter','latex');
hold on
for j = 1:100
    for k=1:21
        [S, Mat, My, E] = rskew5(k-1);
        errors(k) = E;
    end
plot(0:20,0.01*log(errors),'o');
end
%saveas(gca, 'mat508plot.eps','eps');
