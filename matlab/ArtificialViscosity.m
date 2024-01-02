function [U] = ArtificialViscosity(h,U,kappa,e,i,j,nx,ny)


vis         = zeros(nx,ny);
epsilon_p   = zeros(nx,ny);
epsilon_m   = zeros(nx,ny);
epsilon_pp  = zeros(nx,ny);
epsilon_mm  = zeros(nx,ny);

vis(i,j) = abs(h(i+1,j)-2*h(i,j)+h(i-1,j))./...
    (e+abs(h(i+1,j))+2*abs(h(i,j))+abs(h(i-1,j)));


epsilon_p(i,j) = kappa*max(vis(i+1,j),vis(i,j));
epsilon_m(i,j) = kappa*max(vis(i-1,j),vis(i,j));
epsilon_pp(i,j) = kappa*max(vis(i,j+1),vis(i,j));
epsilon_mm(i,j) = kappa*max(vis(i,j-1),vis(i,j));

p(i,j) = U(1,i+1,j)-U(1,i,j);
m(i,j) = U(1,i,j)-U(1,i-1,j);
pp(i,j) = U(1,i,j+1)-U(1,i,j);
mm(i,j) = U(1,i,j)-U(1,i,j-1);

Utemp(i,j) = U(1,i,j);

Utemp(i,j) = Utemp(i,j)+...
    (epsilon_p(i,j).*p(i,j)-...
    epsilon_m(i,j).*(m(i-1,j))+...
    epsilon_pp(i,j).*(pp(i,j))-...
    epsilon_mm(i,j).*(mm(i,j-1)))/2;

U(1,i,j) = Utemp(i,j);