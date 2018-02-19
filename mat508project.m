format long
rng('shuffle')
n=4; k = 10^n;
r =  k*randi([-5,5],1,10);
S = [ 0     r(1) r(2)  r(3)   r(4);...
     -r(1)   0   r(5)  r(6)   r(7);...
     -r(2) -r(5)  0    r(8)   r(9);...
     -r(3) -r(6) -r(8) 0      r(10);...
     -r(4) -r(7) -r(9) -r(10)  0    ]
MatlabExponential = expm(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing my exponential
e = imag(eig(S)); e = e(e>0);
t1 = min(abs(e));
t2 = max(abs(e));

c1 = (t2^3*sin(t1)-t1^3*sin(t2))/(t1*t2*(t2^2-t1^2));
c2 = (t2^4*(1-cos(t1))-t1^4*(1-cos(t2)))/(t1^2*t2^2*(t2^2-t1^2));
c3 = (t2*sin(t1)-t1*sin(t2))/(t1*t2*(t2^2-t1^2));
c4 = (t2^2*(1-cos(t1))-t1^2*(1-cos(t2)))/(t1^2*t2^2*(t2^2-t1^2));

MyExponential = eye(5) + c1*S + c2*S^2 + c3*S^3 + c4*S^4
rel_error = MyExponential-MatlabExponential;
E = norm(rel_error)
