syms x1 x2 x3 x4 a1 a2 a3 b1 b2 c1 c2 c3 c4 d1 d2 r1 r2 s alpha ro N T I M 

x1 = 1.1; x2 = 1.2; x3 = 1.3; x4 = 1.4; 
a1 = 0.05; a2 = 0.15; a3 = 0.1; b1 = 1; b2 = 1; c1 = 1.0; c2 = 0.5; c3 = 1.0; c4 = 1.0; 
d1 = 0.2; d2 = 1.0; r1 = 1.5; r2 = 1.0; s = 0.33; alpha = 0.3; ro = 0.01;

I_dot = s + ro*I*T/(alpha+T) - c1*I*T - d1*I - a1*M*I;

A11 =  -r2*(1+b2*x1);
A12 =  -c4*(x1 + 1/b2);
A13 =  0;
A14 =  -a3*(1/b2 + x1);

A21 =  -c3*x2;
A22 =  r1*(1-b1*x2)-(c2*s/d1 + c3/b2);
A23 =  -c2*x2;
A24 =  -a2*x2;

A31 =  0;
A32 =  -(s + d1*x3)*(alpha*c1 - ro + c1*x2)/(d1*(alpha + x2)); %ro*(x3 + s/d1)/(alpha + x2) -c1*(x3+s/d1) -x4;
% d11 = 0.2; % = d1
A33 =  -d1;%%%%%%%% 0.2
A34 =  -(a1*(s + d1*x3))/d1;

A41 =  0;
A42 =  0;
A43 =  0;
A44 =  -d2;
    
x1_dot = -r2*x1*(1 +b2*x1) -c4/b2*x2 -a3/b2*x4 -c4*x1*x2 -a3*x1*x4;  
x2_dot = r1*x2*(1-b1*x2) -(c2*s/d1 +c3/b2)*x2 -c3*x1*x2 -c2*x2*x3 -a2*x2*x4;
x3_dot = -c2*s/d1*x2 -d1*x3 -a1*s/d1*x4 +(ro*s/d1)*x2/(alpha+x2) +ro*x2*x3/(alpha+x2) -c1*x2*x3 -a1*x3*x4;
x4_dot = -d2*x4 +u;

x33_dot = -d1*x3 -c1*x3*x2 -a1*x3*x4 + ro*x2*(s+d1*x3)/(d1*(alpha+x2)) +(-c1*s*x2 -a1*s*x4)/d1;

% A matrix
A = A31*x1 + A32*x2 + A33*x3 + A34*x4


% x_dot analysis
B = x33_dot
% A,B


A = simplify(A)
B = simplify(B)
% isequaln(A,B)

sub = subs(I_dot,N,x1+1/b2)
sub = subs(sub,T,x2)
sub = subs(sub,I,x3+s/d1)
sub = subs(sub,M,x4)
simplify(sub)

I_to_x3 = -(s + d1*x3)*(alpha*c1 - ro + c1*x2)/(d1*(alpha + x2))*x2 - d1*x3 - ((a1*(s + d1*x3))/d1)*x4;

x3_to_I = s - I*d1 - I*M*a1 - I*T*c1 + (T*c1*s)/d1 - (T*c2*s)/d1 + (I*T*ro)/(T + alpha);

plot(I_to_x3)

% a = [
%    -0.0400   -0.2795   -0.0375;
%    -0.9603   -1.0000    0.0190;
%          0         0   -1.0000]
% az = [
%      0 0        0       0;
%      0 a(1)     a(4)    a(7);
%      0 a(2)     a(5)    a(8);
%      0 a(3)     a(6)    a(9)     ]
