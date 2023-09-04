% Main
prompt ='Escriba el punto que desea ejecutar de la siguiente manera: puntoa, puntob, punto2 o punto3? ';
punto = input(prompt,"s");

n = [5, 10, 20, 50, 100];
dt = [0.1, 0.01, 0.001];

syms x t real;

%%%%%%%%%%%%%%%%%%%%%%%%%% Punto 1a %%%%%%%%%%%%%%%%%%%%%%%%%%
if punto == 'puntoa'
    u_ex(x) = (4/pi)*sin(6*pi*x);
    subplot(1,2,1);
    plot_real(u_ex)
    errores_Linf_1a = zeros(1,length(n));
    errores_L2_1a = zeros(1,length(n));
    for i = 1:length(n)
        [K1a, b1a, a1a, u_h1a, mesh1a, error_L2_1a, error_Linf_1a] = punto1a(n(i), u_ex);
        plot_aprox(u_h1a, mesh1a, n(i))
        errores_L2_1a(1,i) = error_L2_1a;
        errores_Linf_1a(1,i) = error_Linf_1a;
    end
    subplot(1,2,2);
    plot_errores(errores_L2_1a, errores_Linf_1a)
%%%%%%%%%%%%%%%%%%%%%%%%%% Punto 1b %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif punto == 'puntob'
    u_ex = sin(pi*x)*sin(t*pi);
    errores1b = zeros(1,length(dt));
    subplot(1,2,1);
    for i = 1:length(dt)
        [K1b, M1b, b1b, a1b, u_h1b, mesh1b, error1b] = punto1b(10, dt(i), 2, u_ex);
        errores1b(1,i) = error1b;
    end
    subplot(1,2,2);
    plot_errores(errores1b)
%%%%%%%%%%%%%%%%%%%%%%%%%% Punto 2 %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif punto == 'punto2'
    u_ex(x) = x^4 - 2*x^3 + x^2;
    subplot(1,2,1);
    plot_real(u_ex)
    errores_Linf_2 = zeros(1,length(n));
    errores_L2_2 = zeros(1,length(n));
    for i = 1:length(n)
        [K2, M2, P2, Q2, b2, d2, a2, c2, u_h2, mesh2, error_L2_2, error_Linf_2] = punto2(n(i), u_ex);
        plot_aprox(u_h2, mesh2, n(i))
        errores_L2_2(1,i) = error_L2_2;
        errores_Linf_2(1,i) = error_Linf_2;
    end
    subplot(1,2,2);
    plot_errores(errores_L2_2, errores_Linf_2)
%%%%%%%%%%%%%%%%%%%%%%%%%% Punto 3 %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif punto == 'punto3'
    u_ex(x) = x^4 - 2*x^3 + x^2;
    subplot(1,2,1);
    plot_real(u_ex)
    errores_L2_3 = zeros(1,length(n));
    errores_Linf_3 = zeros(1,length(n));
    for i = 1:length(n)
        [K3, M3, P3, Q3, b3, d3, a3, c3, u_h3, mesh3, error_L2_3, error_Linf_3] = punto3(n(i), u_ex);
        plot_aprox(u_h3, mesh3, n(i))
        errores_L2_3(1,i) = error_L2_3;
        errores_Linf_3(1,i) = error_Linf_3;
    end
    subplot(1,2,2);
    plot_errores(errores_L2_3, errores_Linf_3)
end
