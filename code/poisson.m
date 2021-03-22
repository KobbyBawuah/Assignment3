%Assign 3 
%q2
%Kwabena Gyasi Bawuah 
%101048814

% using question 1 a as background 

%initiailizing the dimensions of our matrices, ensuring L is 3/2 times W

W = 100;                
L = (3/2)*W;    

G = sparse(W*L);      
Op = zeros(1, W*L);    

boxD = L*(2/5);
boxU = W*(3/5);
boxR = L*(3/5);
boxL = W*(2/5); 
%Editting Assig 2 part 2a
%sigma reference in and out the box
sigOut = 1; 
sigIn = 10^-2;

% leftEdge = midX - boxL/2;
% rightEdge = midX + boxL/2;
% topEdge = midY + boxW/2;
% bottomEdge = midY - boxW/2;

box = [boxL boxD boxU boxR];


%difine the dimension with the given sigma
for i = 1:W
    for j = 1:L
        if (i > box(1) && i < box(2) && (j < box(3)||j > box(4)))
            sigmamap(i, j) = sigIn;
        else
            sigmamap(i, j) = sigOut;  
        end
    end
end

%then we fill in the whole matrix with the sigma values defined
for i = 1:W 
    for j = 1:L       
        n = j + (i-1)*L;
        nxps = j + (i)*L;
        nxms = j + (i-2)*L;
        nyps = (j + 1) + (i-1)*L;
        nyms = (j - 1) + (i-1)*L;
        
        if (i == 1)
            G(n, :) = 0;
            G(n,n) = 1;
            Op(n) = 1;
            %already assigned above
            %sigmaMap(i,j) = sigOut;
        elseif (i == W)
            G(n, :) = 0;
            G(n,n) = 1;
            Op(n) = 0;
            %sigmaMap(i,j) = sigOut;
        elseif (j == 1)
            G(n,nxps) = (sigmamap(i+1, j) + sigmamap(i,j))/2;
            G(n, nxms) = (sigmamap(i-1, j) + sigmamap(i,j))/2;
            G(n, nyps) = (sigmamap(i, j+1) + sigmamap(i,j))/2;
            G(n,n) = -(G(n,nxps)+G(n,nxms)+G(n,nyps));
        elseif (j == L)
            G(n,nxps) = (sigmamap(i+1, j) + sigmamap(i,j))/2;
            G(n, nxms) = (sigmamap(i-1, j) + sigmamap(i,j))/2;
            G(n, nyms) = (sigmamap(i, j-1) + sigmamap(i,j))/2;
            G(n,n) = -(G(n,nxps)+G(n,nxms)+G(n,nyms));
        else
            G(n,nxps) = (sigmamap(i+1, j) + sigmamap(i,j))/2;
            G(n, nxms) = (sigmamap(i-1, j) + sigmamap(i,j))/2;
            G(n, nyps) = (sigmamap(i, j+1) + sigmamap(i,j))/2;
            G(n, nyms) = (sigmamap(i, j-1) + sigmamap(i,j))/2;
            G(n,n) = -(G(n,nxps)+G(n,nxms)+G(n,nyps)+G(n,nyms));
%             G(n,n) = -3;
%             if(i>leftEdge && i<rightEdge)
%                 G(n,nxms) = sigIn;
%                 G(n,nxps) = sigIn;
%                 G(n,nyms) = sigIn;
%                 sigmaMap(i,j) = sigIn;
%             else
%                 G(n,nxms) = sigOut;
%                 G(n,nxps) = sigOut;
%                 G(n,nyms) = sigOut;
%                 sigmaMap(i,j) = sigOut;
%             end
%         elseif (j == 1)
%             G(n,n) = -3;
%             if(i>leftEdge && i<rightEdge)
%                 G(n,nxms) = sigIn;
%                 G(n,nxps) = sigIn;
%                 G(n,nyps) = sigIn;
%                 sigmaMap(i,j) = sigIn;
%             else
%                 G(n,nxms) = sigOut;
%                 G(n,nxps) = sigOut;
%                 G(n,nyps) = sigOut;
%                 sigmaMap(i,j) = sigOut;
%             end
%         else
%             G(n,n) = -4;
%             if( (j>topEdge || j<bottomEdge) && i>leftEdge && i<rightEdge)
%                 G(n,nxps) = sigIn;
%                 G(n,nxms) = sigIn;
%                 G(n,nyps) = sigIn;
%                 G(n,nyms) = sigIn;
%                 sigmaMap(i,j) = sigIn;
%             else
%                 G(n,nxps) = sigOut;
%                 G(n,nxms) = sigOut;
%                 G(n,nyps) = sigOut;
%                 G(n,nyms) = sigOut;
%                 sigmaMap(i,j) = sigOut;
%             end
        end
    end
end

Voltage = G\Op';
sol = zeros(L, W, 1);

for i = 1:W
    for j = 1:L
        n = j + (i-1)*L;
        sol(j,i) = Voltage(n);
    end
end

%The electric field can be derived from the surface voltage using a
%gradient
[Ey, Ex] = gradient(sol);
%J, the current density, is calculated by multiplying sigma and the
%electric field together. Combing the x and y matrices, a surface plot is
%derived by surfing this matrix.
J_x = sigmamap'.*Ey;
J_y = sigmamap'.*Ex;
J = sqrt(J_x.^2 + J_y.^2);

% Sigma(x,y) Surface Plot
% figure(9)
% subplot(2,1,1);
% surf(sigmamap);
% xlabel('x');
% ylabel('y');
% zlabel('V(x,y)')
% title('Sigma Charge Density Plot');

%subplot(2,1,2);
figure(1)
surf(sol)
xlabel("X position")
ylabel("Y position")
zlabel("Voltage")
title('Voltage Plot');

%X component of electric field surface plot
% figure(3)
% subplot(2,1,1);
% surf(-Ey)
% xlabel('x');
% ylabel('y');
% zlabel('V(x,y)')
% title('Electric Field Plot for x');

%Y component of electric field surface plot
% subplot(2,1,2);
% surf(-Ex)
% xlabel('x');
% ylabel('y');
% zlabel('V(x,y)')
% title('Electric Field Plot for y');

% figure(5)
% surf(J)
% xlabel('x');
% ylabel('y');
% zlabel('V(x,y)')
% title('Current Density σ(x, y)');

figure(2)
quiver(Ex,Ey);
xlabel('x');
ylabel('y');
title("2-D plot")
