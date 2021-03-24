% CEE 279 Project 2 problem 1 Spring 2011
% Jonathan Frame
clear, clc, 
dx          = 25;
dy          = 25;
dt          = 0.1;
tau         = dt/dx/2;
t           = 0;
xmax        = 1000;
ymax        = 500;
tmax        = 1000;
x           = 0:dx:xmax;
y           = 0:dy:ymax;
nx          = length(x);
ny          = length(y);
time        = 0:dt:tmax;
nt          = length(time);
U           = zeros(3,nx,ny);
Us          = zeros(3,nx,ny);
Uss         = zeros(3,nx,ny);
F           = zeros(3,nx,ny);
E           = zeros(3,nx,ny);
lambda      = zeros(3,nx,ny);
S           = zeros(3,nx,ny);
v           = zeros(nx,ny);
h           = zeros(nx,ny);
u           = zeros(nx,ny);
g           = 9.8;
count       = 0.000;
n           = 0.2;     
kappa       = 0.6;
h0          = 5;
e           = 0.00001;
ntplot      = 20; nplot = 0;
x0          = round(nx/2);

h(1:nx,1:ny) = 5;
amp = 2;
U(1,:,:)                = h;
F(3,:,:)                = g*(h(:,:).^2)/2;
E(2,:,:)                = g*(h(:,:).^2)/2;
htot                    = sum(h);
alt = 0;

boundary = h*0;
boundary(round(900/dx):nx,round(100/dy):round(200/dy)) = 1;
boundary(1:round(800/dx),round(100/dy):round(200/dy)) = 1;

boundary(round(800/dx),round(100/dy):round(200/dy)) = 1;
boundary(round(900/dx),round(100/dy):round(200/dy)) = 1;

while t < tmax
    
    t = t+dt; count = count + 1; nplot = nplot+1;
    
    alt = randperm(4)-1;
    alt = alt(1);
    
    wave = h0+amp*sin(t/2/pi);
    
    U(1,:,1)                = wave;
    F(3,:,1)                = g*(wave^2)/2;
    E(2,:,1)                = g*(wave^2)/2;  
    
    i = 2:nx-1;  j = 2:ny-1;
    
    %---------------Predictor---------%
    Us = predictor(U,E,F,S,i,j,alt,tau,dt);
    
    %----------Corrector---------%
    Uss = corrector(Us,E,F,S,i,j,alt,tau,dt);
    
    %-----------True U values---------%
    for i = 2:nx-1
        for j = 2:ny-1
             if boundary(i,j) ~= 1 && boundary(i-1,j) ~= 1 &&...
                    boundary(i+1,j) ~= 1 &&...
                    boundary(i,j+1) ~= 1 && boundary(i,j-1) ~= 1
            U(:,i,j) = (1/2)*(Us(:,i,j)+Uss(:,i,j));
            else
                if i>2 && j > 2 && i < nx-1 && j < ny-1
                    if boundary(i,j-1) ~= 1 && boundary(i-1,j) == 1 &&...
                    boundary(i+1,j) == 1 && boundary(i,j+1) == 1 
                        U(2:3,i,j-1) = -U(2:3,i,j-2);
                        U(1,i,j-1) = U(1,i,j-2);
                    end 
                    if boundary(i,j+1) ~= 1 && boundary(i-1,j) == 1 &&...
                    boundary(i+1,j) == 1 && boundary(i,j-1) == 1
                        U(2:3,i,j+1) = -U(2:3,i,j+2);
                        U(1,i,j+1) = U(1,i,j+2);
                    end
                    if boundary(i-1,j) ~= 1 && boundary(i+1,j) == 1 &&...
                    boundary(i,j+1) == 1 && boundary(i,j-1) == 1
                        U(2:3,i-1,j) = -U(2:3,i-2,j);
                        U(1,i-1,j) = U(1,i-2,j);
                    end 
                    if boundary(i+1,j) ~= 1 && boundary(i-1,j) == 1 &&...
                    boundary(i,j+1) == 1 && boundary(i,j-1) == 1
                        U(2:3,i+1,j) = -U(2:3,i+2,j);
                        U(1,i+1,j) = U(1,i+2,j);  
                    end                
                end
            end
            if boundary(i,j) == 1 && boundary(i-1,j) == 1 &&...
                    boundary(i+1,j) == 1 &&...
                    boundary(i,j+1) == 1 && boundary(i,j-1) == 1
                U(:,i,j) = h0+1;
            end
        end
    end
    
    %-----------Artificial Viscosity---------------%
    U = ArtificialViscosity(h,U,kappa,e,i,j,nx,ny);
    
    %------------ Boundary Conditions---------------%
    U(1,:,:)    = max(0,U(1,:,:));
    U(1,:,1)  = U(1,:,2);
    U(1,nx,:)   = U(1,nx-1,:);  U(1,:,ny) = U(1,:,ny-1);
    U(2,:,1)  = -U(2,:,2);
    U(2,nx,:)   = -U(2,nx-1,:);  U(2,:,ny) = -U(2,:,ny-1);
    U(3,:,1)  = -U(2,:,2);
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
        mass = sum(h)/htot,t
        figure(1)
        pcolor(x,y,h'), 
        xlabel('x')
        ylabel('y')
%        zlabel('water height')
        figure(2)
        plot(h(round(900/dx),:))
        axis([0 ny -amp+h0 amp+h0])  
        pause(0.001)
        nplot = 0;
    end

end
