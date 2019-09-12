f = importdata('C:\Users\kosta\Desktop\Cholesterol2.dat');
% scatter(f(:,1),f(:,2))
Chol = f(:,2); %Cholesterol decrease
Comp = f(:,1); %Compliance

A_cub = [ones(size(Chol)),Comp,Comp.^2,Comp.^3];
A_quad = [ones(size(Chol)),Comp,Comp.^2,Comp.^3,Comp.^4];

x_cub = A_cub\Chol; %lsqr(A_cub,Chol);
x_quad = A_quad\Chol; %lsqr(A_quad,Chol);

% f_c = @(t) x_cub(1)+x_cub(2)*t+x_cub(3)*t^2+x_cub(4)*t^3;
% f_q = @(t) x_quad(1)+x_quad(2)*t+x_quad(3)*t^2+x_quad(4)*t^3+x_quad(5)*t^4;
%======================== Cubic Model =========================================
t = linspace(-2.5,2.5,1000);
model_cub = x_cub(1)+x_cub(2)*t+x_cub(3)*t.^2+x_cub(4)*t.^3;
% figure(1);
% plot(t,model_cub)
% hold on
% scatter(Comp,Chol)
% title('a_0+a_1*x+a_2*x^2+a_3*x^3');
% xlabel('Compliance level');
% ylabel('Cholesterol value');
% hold off
%======================== Quartic Model =========================================
model_quad = x_quad(1)+x_quad(2)*t+x_quad(3)*t.^2+x_quad(4)*t.^3+x_quad(5)*t.^4;
% figure(2);
% plot(t,model_quad)
% hold on
% scatter(Comp,Chol)
% title('a_0+a_1*x+a_2*x^2+a_3*x^3+a_4*x^4');
% xlabel('Compliance level');
% ylabel('Cholesterol value');
% hold off
%================================================================================
Err(1) = norm(Chol - (A_cub*x_cub));
Err(2) = norm(Chol - (A_quad*x_quad));
disp(Err)



%========================= TLS ========================================
Atls = [A_cub,Chol];
[U,S,V] = svd(Atls,'econ');
n = size(Atls,2)-1;
X = -(1/V(n+1,n+1))*V(:,n+1);
Xtls = X(1:end-1);

mod_TLS = Xtls(1)+Xtls(2)*t+Xtls(3)*t.^2+Xtls(4)*t.^3;
% figure(2)
% plot(t,mod_TLS)
% hold on
% scatter(Comp,Chol)

% Delta = -Atls*V(:,n+1)*V(:,n+1)';
Delta = -Atls*V(:,end)*V(:,end)';
% Delta2 = -S(n+1,n+1)*U(:,n+1)*V(:,n+1)';
Delta_TLS = norm(Delta,'fro');
% Delta2_LST = norm(Delta2,'fro');
%======================= PLOTS ===========================================
figure(1)
hold on
plot(t,model_cub,'g','DisplayName','Cubic Model')
plot(t,model_quad,'b','DisplayName','Quartic Model')
plot(t,mod_TLS,'r','DisplayName','TLS')
scatter(Comp,Chol,'k','DisplayName','Data')
xlabel('Compliance level');
ylabel('Cholesterol value');
title('Fitting Models');
% legend({'a_0+a_1*x+a_2*x^2+a_3*x^3','a_0+a_1*x+a_2*x^2+a_3*x^3+a_4*x^4','Total Least Squares','Scatter'},'Location','southwest')
legend('location','southwest')
hold off
%=========================================================================
if (S(n,n)>S(n+1,n+1))
    disp('Monadki lysh!')
else
    disp('Maybe there is no solution!')
end
