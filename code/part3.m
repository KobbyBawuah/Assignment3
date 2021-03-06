%Assign 3 
%q3
%editting coupled
close all;
clc
%Kwabena Gyasi Bawuah
%101048814
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%electron spec
 global C

    addpath ../geom2d/geom2d

    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                    % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665; %metres (32.1740 ft) per s²
    
    T = 300;
    k = 1.38e-23;
    mn = 0.26*C.m_0; %effective mass
    tmn = 0.2e-12;    % Mean time between collisions
   
    vth = sqrt((2*C.kb*T)/mn);% Thermal velocity
    
    freepath = vth*tmn   % mean free path
    
    ConductorL = 180e-9;
    ConductorW = 80e-9;
    
    dpoints = 5e4;
    ecount = 15; %the number of electron to show on plot
    
    %electron concentratoin 
    den = 1e15*100^-2;
    
    detaT= ConductorW/vth/100;
    sims = 200;
    
%     Xpos = rand(1,ecount).*ConductorW;
%     Ypos = rand(1,ecount).*ConductorL;
    traj=zeros(sims,ecount*2);
    temp=zeros(sims, 1);
    
    temp(:,1)= 300;
    
    
    MFP = vth*0.2e-12;
    
    Vx = 0.8;
    Vy = 0;
    dens = 1e15*100^-2;
    
    Ex = Vx/ConductorL
    Ey = Vy/ConductorW
    
    Fx = -C.q_0*Ex
    Fy = -C.q_0*Ey
    
    dVx = Fx*detaT/mn;
    dVy = Fy*detaT/mn;
    dVx = dVx.*ones(dpoints,1);
    dVy = dVy.*ones(dpoints,1);
    
    Pscat = 1-exp(-detaT/tmn);
    %----------------------------------------
    ProbDistr = makedist('Normal','mu', 0, 'sigma', sqrt(C.kb*T/mn));
    
    tspec = 0;
    bspec=0;
    
    
    boxes = 1e-9.*[80 120 0 40; 80 120 60 100];
    specularbox = [0 1];

    for i = 1: dpoints
        angle = rand*2*3.14;
        state(i,:)= [ConductorL*rand ConductorW*rand random(ProbDistr) random(ProbDistr)];
        
        if (state(i,2)>60e-9 &(state(i,1)>80e-9 & state(i,1)<120e-9))  | (state(i,2)< 40e-9 &(state(i,1)>80e-9 & state(i,1)<120e-9))
            state(i,1:2) = [ConductorL*rand ConductorW*rand];
        end
 
    end
    
%     %take plot of loop to stop the refresh
% figure(5);
% subplot(4,1,2);
% %         plot(detaT*(0:i-1), temp(1:i));
% %         plot(detaT*(0:i-1),temp(1:i));
% tPlot = animatedline;
% xlabel('time(s)');
% ylabel('Temperature (K)');
% title('Temperature of semiconductor over time');
% 
% figure(5);        
% subplot(4,1,3);
% %part = sqrt(state(:,3).^2 + state(:,4).^2);
% % xlim([0 7e5]);
% % ylim([0 2000]);
% %histogram(part);
% currentPlot = animatedline ;
% xlabel('t(s)');
% ylabel('Current density');
% title('Drift');
    
    for i = 1:sims
        
        state(:,3) = state(:,3) + dVx;
        state(:,4) = state(:,4) + dVy;
        state(:,1:2) = state(:,1:2) + detaT.*state(:,3:4);
    
        out = state(:,1) > ConductorL;
        state(out,1) = state(out,1) - ConductorL;

        out = state(:,1) < 0;
        state(out,1) = state(out,1) + ConductorL;

        out = state(:,2) > ConductorW;

    if(tspec)
        state(out,2) = 2*ConductorW - state(out,2);
        state(out,4) = -state(out,4);
    else 
        state(out,2) = ConductorW;
        part = sqrt(state(out,3).^2 + state(out,4).^2);
        angle = rand([sum(out),1])*2*3.14;
        state(out,3) = part.*cos(angle);
        state(out,4) = -abs(part.*sin(angle));
    end
    
    out = state(:,2) < 0;
    
    if(bspec)
        state(out,2) = -state(out,2);
        state(out,4) = -state(out,4);
    else 
        state(out,2) = 0;
        part = sqrt(state(out,3).^2 + state(out,4).^2);
        angle = rand([sum(out),1])*2*3.41;
        state(out,3) = part.*cos(angle);
        state(out,4) = abs(part.*sin(angle));
    end
    for out=1: dpoints
            if (state(out,2)>60e-9 &(state(out,1)>80e-9 & state(out,1)<120e-9)) 
                boxNum = 1;
            elseif (state(out,2)< 40e-9 &(state(out,1)>80e-9 & state(out,1)<120e-9))
                boxNum = 2;
            else 
                boxNum = 0;
            end
            while(boxNum ~= 0)
                XDist = 0;
                newx = 0;
                if(state(out,3) > 0)
                    XDist = state(out,1) - boxes(boxNum,1);
                    newx = boxes(boxNum,1);
                else
                    XDist = boxes(boxNum,2) - state(out,1);
                    newx = boxes(boxNum,2);
                end
                yDist = 0;
                newy = 0;
                if(state(out,4) > 0)
                    yDist = state(out,2) - boxes(boxNum, 3);
                    newy = boxes(boxNum, 3);
                else
                    yDist = boxes(boxNum, 4) - state(out,2);
                    newy = boxes(boxNum, 4);
                end

                if(XDist < yDist)
                    state(out,1) = newx;
                    if(~specularbox(boxNum))
                        sgn = -sign(state(out,3));
                        part = sqrt(state(out,3).^2 + state(out,4).^2);
                        angle = rand()*2*3.14;
                        state(out,3) = sgn.*abs(part.*cos(angle));
                        state(out,4) = part.*sin(angle);
                    else 
                        state(out,3) = -state(out,3);
                    end
                else
                    state(out,2) = newy;
                    if(~specularbox(boxNum))
                        sgn = -sign(state(out,4));
                        part = sqrt(state(out,3).^2 + state(out,4).^2);
                        angle = rand()*2*3.14;
                        state(out,3) = part.*cos(angle);
                        state(out,4) = sgn.*abs(part.*sin(angle));
                    else
                        state(out,4) = -state(out,4);
                    end
                end
                 boxNum = 0;
            end
     end

        out = rand(dpoints, 1) < Pscat;
        state(out,3:4) = random(ProbDistr, [sum(out),2]);

        temp(i) = (sum(state(:,3).^2) + sum(state(:,4).^2))*mn/k/2/dpoints;
    

        for out=1:ecount
            traj(i, (2*out):(2*out+1)) = state(out,1:2);
        end 
        
        %getting j using density
        %add the points
        
        if  mod(i,5) == 0
        figure(5);
        hold off;
        plot(state(1:ecount,1)./1e-9, state(1:ecount,2)./1e-9, 'o');
        hold on;
        
        for ouut=1:size(boxes,1)
           plot([boxes(ouut, 1) boxes(ouut, 1) boxes(ouut, 2) boxes(ouut, 2) boxes(ouut, 1)]./1e-9,[boxes(ouut, 3) boxes(ouut, 4) boxes(ouut, 4) boxes(ouut, 3) boxes(ouut, 3)]./1e-9, 'k-');
        end
        
        xlim([0 ConductorL/1e-9])
        ylim([0 ConductorW/1e-9])
        xlabel('Electrons x position (nm)')
        ylabel('Electrons y position (nm)')
        title('Electrons Simulations')
        
%         figure(6)
%         part = sqrt(state(:,3).^2 + state(:,4).^2);
%         xlim([0 7e5]);
%         ylim([0 2000]);
%         histogram(part);
%         xlabel('v(m/s)');
%         ylabel('Particle count');
%         title('Histogram to show particle speed');
        
        J(i, 1) = -C.q_0.*dens.*mean(state(:,3));
        J(i, 2) = -C.q_0.*dens.*mean(state(:,4))
        
        %not needed
%         addpoints(tPlot, detaT.*i, temp(i));
%         addpoints(currentPlot, detaT.*i, J(i,1));
        end
    end 
    
    figure(2)
    hold on;
    xlim([0 ConductorL/1e-9]);
    ylim([0 ConductorW/1e-9]);
    xlabel('Electrons x position (nm)');
    ylabel('Electrons y position (nm)');
    title('Trajectories of Electrons with curl from field');

    for i = 1: ecount
        plot(traj(:,i*2)./1e-9, traj(:,i*2+1)./1e-9, '-');
    end
   
   for out=1:size(boxes,1)
   plot([boxes(out, 1) boxes(out, 1) boxes(out, 2) boxes(out, 2) boxes(out, 1)]./1e-9,...
       [boxes(out, 3) boxes(out, 4) boxes(out, 4) boxes(out, 3) boxes(out, 3)]./1e-9, 'k-');
    end

%     figure(7)
%     part = sqrt(state(:,3).^2 + state(:,4).^2);
%     xlim([0 7e5]);
%     ylim([0 2000]);
%     histogram(part);
%     xlabel('v(m/s)');
%     ylabel('Particle count');
%     title('Histogram to show particle speed');

    dens = hist3(state(:,1:2),[200 100])';

    N = 20;
    sigma = 3;
    [x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
    f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
    f=f./sum(f(:));
    figure(3);
    imagesc(conv2(dens,f,'same'));
    set(gca,'YDir','normal');
    title('Electron Density');
    xlabel('Electrons x position (nm)');
    ylabel('Electrons y position (nm)');

%     tempSumX = zeros(ceil(ConductorL/1e-9),ceil(ConductorW/1e-9));
%     tempSumY = zeros(ceil(ConductorL/1e-9),ceil(ConductorW/1e-9));
%     tempSum = zeros(ceil(ConductorL/1e-9),ceil(ConductorW/1e-9));
% 
%     for i=1:dpoints
% 
%         x = floor(state(i,1)/1e-9);
%         y = floor(state(i,2)/1e-9);
%         if(x==0)
%             x = 1;
%         end
%         if(y==0)
%             y= 1;
%         end
% 
%         tempSumY(x,y) = tempSumY(x,y) + state(i,3)^2;
%         tempSumX(x,y) = tempSumX(x,y) + state(i,4)^2;
%         tempSum(x,y) = tempSum(x,y) + 1;
%     end
% 
%     temp = (tempSumX + tempSumY).*mn./k./2./tempSum;
%     temp(isnan(temp)) = 0;
%     temp = temp';
% 
%     N = 20;
%     sigma = 3;
%     [x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
%     f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
%     f=f./sum(f(:));
%     figure(9);
%     surf(conv2(temp,f,'same'));
%     %set(gca,'YDir','normal');
%     title('Map of temperature');
%     xlabel('Electrons x position (nm)');
%     ylabel('Electrons y position (nm)');
    