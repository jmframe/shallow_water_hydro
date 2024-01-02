% CEE 279 Project 2 problem 1 Spring 2011
% Jonathan Frame
%clear, clc, 
dx          = 1;
dy          = 1;
dt          = 0.001;
tau         = dt/dx/2;
t           = 0;
xmax        = 50;
ymax        = 20;
tmax        = 30;

[x,y,nx,ny,time,nt,U,Us,Uss,F,E,lambda,S,v,h,u] =...
    make_vectors(dx,dy,dt,xmax,ymax,tmax);

g           = 9.8;
kappa       = 0.6;
count       = 0.000;
n           = 0.00;     
h0          = 4;
h1          = 2;
e           = 0.00001;
ntplot      = 15; nplot = 0;
x0          = round(nx/2);

for i = 1:nx
    for j = 1:ny
        h(i,j) = h1*exp(-(x(i)^2+y(j)^2)/4)+h0;
    end
end

U(1,:,:)                = h;
F(3,:,:)                = g*(h(:,:).^2)/2;
E(2,:,:)                = g*(h(:,:).^2)/2;
htot                    = sum(h);
alt = 0;

while t < tmax
    
    t = t+dt; count = count + 1; nplot = nplot+1;
    
    alt = randperm(4)-1;
    alt = alt(1);
    
    i = 2:nx-1;  j = 2:ny-1;
    
    %---------------Predictor---------%
    Us = predictor(U,E,F,S,i,j,alt,tau,dt);
    
    %----------Corrector---------%
    Uss = corrector(Us,E,F,S,i,j,alt,tau,dt);
    
    %-----------True U values---------%
    U(:,i,j) = (1/2)*(Us(:,i,j)+Uss(:,i,j));
    
    %-----------Artificial Viscosity---------------%
    U = ArtificialViscosity(h,U,kappa,e,i,j,nx,ny);
    
    %------------ Boundary Conditions---------------%
    U(1,:,:)    = max(0,U(1,:,:));
    U(1,1,:)    = U(1,2,:);     U(1,:,1)  = U(1,:,2);
    U(1,nx,:)   = U(1,nx-1,:);  U(1,:,ny) = U(1,:,ny-1);
    U(2,1,:)    = -U(2,2,:);     U(2,:,1)  = -U(2,:,2);
    U(2,nx,:)   = -U(2,nx-1,:);  U(2,:,ny) = -U(2,:,ny-1);
    U(3,1,:)    = -U(3,2,:);     U(3,:,1)  = -U(2,:,2);
    U(3,nx,:)   = -U(3,nx-1,:);  U(3,:,ny) = -U(2,:,ny-1);

    i = 1:nx;  j = 1:ny;
    
    %----------Solving for h and u,v and F,E and S-------%
    h(i,j) = U(1,i,j);
    u(i,j) = U(2,i,j)./(U(1,i,j)+e);
    v(i,j) = U(3,i,j)./(U(1,i,j)+e);
        
    E(1,:,:) = u.*h;    
    E(2,:,:) = h.*u.^2+g*h.^2/2;  
    E(3,:,:) = u.*v.*h;
    F(1,:,:) = v.*h;    
    F(2,:,:) = u.*v.*h;  
    F(3,:,:) = h.*v.^2+g*h.^2/2;
    
    S(2,:,:) = g*h.*(abs(u).*u.*n^2./(e+h.^(1/3)));
    S(3,:,:) = g*h.*(abs(v).*v.*n^2./(e+h.^(1/3)));

    if nplot == ntplot
        mass = sum(h)/htot,t,
        surf(y,x,h), 
       axis([-ymax ymax -xmax xmax h0-h1 h0+h1])
        xlabel('x')
        ylabel('y')
       zlabel('water height')
        pause(0.001)
        nplot = 0;
    end

end
