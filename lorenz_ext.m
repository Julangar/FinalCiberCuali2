function f=lorenz_ext(t,X)
%
%  Lorenz equation 
%
%               dx/dt = SIGMA*(y - x)
%               dy/dt = R*x - y -x*z
%               dz/dt= x*y - BETA*z
%
%        In demo run SIGMA = 10, R = 28, BETA = 8/3
%        Initial conditions: x(0) = 0, y(0) = 1, z(0) = 0;
%        Reference values for t=10 000 : 
%              L_1 = 0.9022, L_2 = 0.0003, LE3 = -14.5691
%
%        See:
%    K. Ramasubramanian, M.S. Sriram, "A comparative study of computation 
%    of Lyapunov spectra with different algorithms", Physica D 139 (2000) 72-86.
%
% --------------------------------------------------------------------
% Copyright (C) 2004, Govorukhin V.N.


% Values of parameters
a = 0.4;
v = 0.5;

x=X(1); y=X(2); z=X(3);

Y= [X(4), X(7), X(10);
    X(5), X(8), X(11);
    X(6), X(9), X(12)];

f=zeros(9,1);

%Lorenz equation
f(1)=-v*x+2*y;
f(2)=-v*y+(z-a)*x;
f(3)=1-x*y;


%Linearized system

 Jac=[-v, 2,     0;
         z-a,    -v,    x;
           y,     x, 0];
  
%Variational equation   
f(4:12)=Jac*Y;

%Output data must be a column vector


