function [K, M, b, a, u_h, mesh_espacial, error_L2] = punto1b(n, dt, T, u_ex)

syms x t real
% Definir funcion
dx1u_ex = diff(u_ex, x);
dx2u_ex = diff(dx1u_ex, x);
dt1u_ex = diff(u_ex, t);
f = dt1u_ex-dx2u_ex+u_ex;

% Dominio espacial
x_inf = 0;
x_sup = 1;
% Dominio temporal
t_inf = 0; 
t_sup = T;

N = 2*n-1; 
hx = 1/n; % tamaño de paso espacial
ht = dt; % tamaño de paso temporal
Nt = T/ht;
mesh_espacial = x_inf:hx/2:x_sup;

% Inicialización de matriz K y vector b
K = zeros(N, N);
M = zeros(N, N);
a = zeros(N,Nt);

% Funciones phi y psi
syms w z positive integer

phi(x) = piecewise((-1 <= x) & (x <= 0), (1+x)*(1+2*x), (0 < x) & (x <= 1), (1-x)*(1-2*x), abs(x) > 1, 0);
psi(x) = piecewise(abs(x) <= 1/2, 1-4*x.^2, abs(x) > 1/2, 0);

y = (x-w.*(hx/2))/(hx);

% Función de prueba lambda
lambda(z) = piecewise(mod(z,2) == 0, phi(subs(y,w,z)), mod(z,2) == 1, psi(subs(y,w,z)));

nodos = 1:N;
v = lambda(nodos);
dv = diff(v);

% Se llenan las matrices
for i = 1:N
    for j=i:min(i+2,N)
       K(i,j) = eval(int(dv(i)*dv(j),x_inf,x_sup));
       K(j,i) = K(i,j);
       M(i,j) = eval(int(v(i)*v(j),x_inf,x_sup));
       M(j,i) = M(i,j);
    end
end

% Calcular el vector a(t=0) y el vector b
u_exacta(t,x) = u_ex;
u_0(x) = u_exacta(0,x);
prod_pto = eval(int(u_0*v,x_inf,x_sup))';
a(:,1) = M\prod_pto;

b = eval(int(f*v,x_inf,x_sup))';
B(t) = b;
P = M + ht*K + ht*M;
for k=1:Nt
    Q = ht*B(k*ht) + M*a(:,k);
    a(:,k+1) = P\Q;
end

% Se obtiene la u aproximada
u_h = v*a;

% Se obtiene el error L2
errores_L2 = zeros(Nt,1);
mesh_temporal = t_inf:ht:t_sup;
for k=1:Nt
    m = mesh_temporal(k);
    errores_L2(k) = double(int((u_h(k)-u_exacta(m, x))^2,x,[0,1])^1/2);
end
error_L2 = mean(errores_L2);

time_plot(u_exacta, u_h, ht, mesh_espacial)

end