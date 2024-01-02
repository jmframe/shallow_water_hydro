function [Uss] = corrector(Us,E,F,S,i,j,alt, tau,dt)
    if alt == 0
        Uss(:,i,j) = Us(:,i,j)-...
            tau*(E(:,i+1,j)-E(:,i,j))-tau*(F(:,i,j+1)-F(:,i,j))-dt*S(:,i,j);
    elseif alt ==1
        Uss(:,i,j) = Us(:,i,j)-...
            tau*(E(:,i,j)-E(:,i-1,j))-tau*(F(:,i,j)-F(:,i,j-1))-dt*S(:,i,j);
    elseif alt ==2
        Uss(:,i,j) = Us(:,i,j)-...
            tau*(E(:,i+1,j)-E(:,i,j))-tau*(F(:,i,j)-F(:,i,j-1))-dt*S(:,i,j);    
    elseif alt ==3
        Uss(:,i,j) = Us(:,i,j)-...
            tau*(E(:,i,j)-E(:,i-1,j))-tau*(F(:,i,j+1)-F(:,i,j))-dt*S(:,i,j);        
    end