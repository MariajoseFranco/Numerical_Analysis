function time_plot(u_ex, u_h, ht, mesh)

syms x t real;
% Graficar solución exacta vs aproximada en distintos valores de t
t_step = 1;
j=1;
for T=0:ht:1
    subplot(2,3,j);
    fplot(u_ex(T, x), [0 1], 'DisplayName', sprintf('u_{ex}(t=%.2f)', T))
    hold all
    U_h(t,x) = u_h(round(t_step,0));
    plot(mesh, eval(subs(U_h(T, x),mesh)), 'DisplayName', sprintf('u_h(t=%.2f)', T))
    title(sprintf('Exacta vs Aproximada (n=%.2f)', 5));
    t_step = t_step + 1;
    j = j + 1;
end

end