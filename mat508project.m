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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%determinant evaluation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
N=100;
dets = zeros(1,N);
Dets = zeros(1,N);

for j = 1:N
    r=randi([-10 10],1,10);
    [M m Er] = expskew5(r,0);
    dets(j)=abs(1-det(m));
    Dets(j)=abs(1-det(M));
end
plot(1:N, dets, 'b',1:N, Dets,'r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% expskew5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MyExponential, MatlabExponential, rel_error] = expskew5(r,k)
% r is a 1x10 vector, the terms correspond to the entries above the
% diagonal of the skew-symmetric matrix S.
format long
S = [ 0     r(1) r(2)  r(3)   r(4);...
     -r(1)   0   r(5)  r(6)   r(7);...
     -r(2) -r(5)  0    r(8)   r(9);...
     -r(3) -r(6) -r(8) 0      r(10);...
     -r(4) -r(7) -r(9) -r(10)  0    ];
k=10^k;
MatlabExponential = expm(k*S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing my exponential
A = charpoly(S); a = A(3); b = A(5);
t1 = sqrt((a-sqrt(a^2-4*b))/2);
t2 = sqrt((a+sqrt(a^2-4*b))/2);

c1 = (t2^3*sin(k*t1)-t1^3*sin(k*t2))/(t1*t2*(t2^2-t1^2));
c2 = (t2^4*(1-cos(k*t1))-t1^4*(1-cos(k*t2)))/(t1^2*t2^2*(t2^2-t1^2));
c3 = (t2*sin(k*t1)-t1*sin(k*t2))/(t1*t2*(t2^2-t1^2));
c4 = (t2^2*(1-cos(k*t1))-t1^2*(1-cos(k*t2)))/(t1^2*t2^2*(t2^2-t1^2));

MyExponential = eye(5) + c1*S + c2*S^2 + c3*S^3 + c4*S^4;
E = MyExponential-MatlabExponential;
rel_error = norm(E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%presentation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [interestingmatrix, I, notsointerestingmatrix, i] = determinantanalysis(N,s) 
Matrices = zeros(5,5*N);
close all
dets = zeros(1,N);
Dets = zeros(1,N);

for j = 1:N
    r=randi([-s s],1,10);
    [M, m, ~, S] = expskew5(r,0);
    dets(j)=abs(1-det(m));
    Dets(j)=abs(1-det(M));
    Matrices(:,5*j-4:5*j) = S;
end
[~, I] = max(abs(dets-Dets)); I=I(1);
[~, i] = min(abs(dets-Dets)); i=i(1);
interestingmatrix = Matrices(:,5*I-4:5*I);
notsointerestingmatrix = Matrices(:,5*i-4:5*i);
plot(1:N, dets, '--ob',1:N, Dets,'-*r');
xlabel('Generated Matrices S')
ylabel('$|1-\det(S)|$','Interpreter','latex')
legend('Matlab''s expm','Closed Formula')


