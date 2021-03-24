function [Us] = predictor(U,E,F,S,i,j,alt, tau,dt)
    if alt ==0
        Us(:,i,j) = U(:,i,j)-...
            tau*(E(:,i+1,j)-E(:,i,j))-tau*(F(:,i,j)-F(:,i,j-1))-dt*S(:,i,j); 
    elseif alt ==1
        Us(:,i,j) = U(:,i,j)-...
            tau*(E(:,i,j)-E(:,i-1,j))-tau*(F(:,i,j+1)-F(:,i,j))-dt*S(:,i,j);   
    elseif alt == 2
        Us(:,i,j) = U(:,i,j)-...
            tau*(E(:,i+1,j)-E(:,i,j))-tau*(F(:,i,j+1)-F(:,i,j))-dt*S(:,i,j);
    elseif alt ==3
        Us(:,i,j) = U(:,i,j)-...
            tau*(E(:,i,j)-E(:,i-1,j))-tau*(F(:,i,j)-F(:,i,j-1))-dt*S(:,i,j);   
    end
