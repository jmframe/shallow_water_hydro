function [x,y,nx,ny,time,nt,U,Us,Uss,F,E,lambda,S,v,h,u] =...
    make_vectors(dx,dy,dt,xmax,ymax,tmax)
    
x           = -xmax:dx:xmax;
y           = -ymax:dy:ymax;
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