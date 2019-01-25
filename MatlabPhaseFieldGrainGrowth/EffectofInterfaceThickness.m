%Captures the Effect of Interface Thickness

clc; clear; close all;

subfldr={'/GG_Q2/'};
sf_dr=subfldr{1};

dt=0.01;
gs=1; %grid size
timesteps=5e6;
plotMOD=25;
MG=200; %number of cells to divide width/height

GRAINS=2;

middle_grain=floor(MG/2);

[xm, ym]=meshgrid(-MG/2:gs:MG/2); %phi locations

[cI, rI] = meshgrid(1:MG+1, 1:MG+1);
circleLogical1 = (rI - middle_grain).^2 + (cI - middle_grain).^2 <= 5.5.^2;

phi(:,:,1)=xm.*0;
phi(:,:,2)=xm.*0+1;

for i=1:MG+1
    for j=1:MG+1
      
      if circleLogical1(i,j)
          phi(i,j,1)=1;
          phi(i,j,2)=0;
      end
      
    end
end

philong=zeros((MG+1)^2, GRAINS);

for g=1:GRAINS
    philong(:,g)=flipud(reshape(flipud(phi(:,:,g)'),(MG+1)^2,1));
end


save([ pwd sf_dr '/InputInformation.mat']);

%phi(:,:,2)-flipud(reshape(philong(:,2)', MG+1, MG+1)')

nnlong=[1:(MG+1)^2]';
NN=flipud(reshape(nnlong', MG+1, MG+1)'); %node numbering

Nx=MG+1;

LapOp=spalloc(Nx*Nx-Nx,Nx*Nx-Nx,5*(Nx-3)*Nx+Nx);

for j=1:length(ym)
 for i=length(xm):-1:1    
     right=j+1;
     if right>size(xm,1)
         right=1;
     end     
     left=j-1;
     if left==0
         left=size(xm,1);
     end     
     down=i+1;
     if down>size(ym,1)
         down=1;
     end     
     up=i-1;
     if up==0
         up=size(ym,1);
     end          
     LapOp(NN(i,j),NN(i,j))=-4;
     LapOp(NN(i,j),NN(i,right))=1;
     LapOp(NN(i,j),NN(i,left))=1;
     LapOp(NN(i,j),NN(down,j))=1;
     LapOp(NN(i,j),NN(up,j))=1;
 end
end
LapOp=LapOp./gs^2;

scount=0;
ff1=figure('Visible','off');
contourf(xm,ym,flipud(reshape(philong(:,1)', MG+1, MG+1)'), 'linestyle','none'); colorbar;
caxis([0 1])
title(['Time is: ' num2str(0)],'Interpreter','latex','FontSize',16)
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
axis equal;
saveas(ff1,[pwd sf_dr sprintf('FIG%d_G1.png',scount)]);
savefig(ff1,[pwd sf_dr sprintf('FIG%d_G1.fig',scount)]);
ff2=figure('Visible','off');
contourf(xm,ym,flipud(reshape(philong(:,2)', MG+1, MG+1)'), 'linestyle','none'); colorbar;
caxis([0 1])
title(['Time is: ' num2str(0)],'Interpreter','latex','FontSize',16)
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
axis equal;
saveas(ff2,[pwd sf_dr sprintf('FIG%d_G2.png',scount)]);
savefig(ff2,[pwd sf_dr sprintf('FIG%d_G2.fig',scount)]);
scount=scount+1;

othergrainindex=[2 1];

middlerow=philong((MG+1)*floor(MG/2):(MG+1)*floor(MG/2)+MG+1,1);

startCircleInd=0;
endCircleInd=0;

for i=2:length(middlerow)-1
    if middlerow(i)<=0.5 && middlerow(i+1)>=0.5
        startCircleInd=i;
    end
    
    if middlerow(i)>=0.5 && middlerow(i+1)<=0.5
        endCircleInd=i;
    end
end
startCircleValue=(0.5-middlerow(startCircleInd))/(middlerow(startCircleInd+1)-middlerow(startCircleInd))+startCircleInd;
endCircleValue=(0.5-middlerow(endCircleInd))/(middlerow(endCircleInd+1)-middlerow(endCircleInd))+endCircleInd;

disp(['phi_1 circle size is ' num2str(endCircleValue-startCircleValue)])

aa=[0 endCircleValue-startCircleValue];


for k=1:timesteps
    
    for g=1:GRAINS
        dfdphi=-philong(:,g)+philong(:,g).^3+3*philong(:,g).*philong(:,othergrainindex(g)).^2;
        philong(:,g)=philong(:,g)-(dt).*(dfdphi-(2).*(LapOp*philong(:,g)));
    end
    
    if mod(k,plotMOD)==0
        ff1=figure('Visible','off');
        contourf(xm,ym,flipud(reshape(philong(:,1)', MG+1, MG+1)'), 'linestyle','none'); colorbar;
        caxis([0 1])
        title(['Time is: ' num2str(k*dt)],'Interpreter','latex','FontSize',16)
        xlabel('X','Interpreter','latex','FontSize',16)
        ylabel('Y','Interpreter','latex','FontSize',16)
        axis equal;
        %saveas(ff1,[pwd sf_dr sprintf('FIG%d_G1.png',scount)]);
        %savefig(ff1,[pwd sf_dr sprintf('FIG%d_G1.fig',scount)]);
        ff2=figure('Visible','off');
        contourf(xm,ym,flipud(reshape(philong(:,2)', MG+1, MG+1)'), 'linestyle','none'); colorbar;
        caxis([0 1])
        title(['Time is: ' num2str(k*dt)],'Interpreter','latex','FontSize',16)
        xlabel('X','Interpreter','latex','FontSize',16)
        ylabel('Y','Interpreter','latex','FontSize',16)
        axis equal;
        %saveas(ff2,[pwd sf_dr sprintf('FIG%d_G2.png',scount)]);
        %savefig(ff2,[pwd sf_dr sprintf('FIG%d_G2.fig',scount)]);
        %save([ pwd sf_dr '/PhiStore' num2str(scount) '.mat'], 'philong','-v7.3');
        scount=scount+1;  
        
        
        middlerow=philong((MG+1)*floor(MG/2):(MG+1)*floor(MG/2)+MG+1,1);

        startCircleInd=0;
        endCircleInd=0;

        for i=2:length(middlerow)-1
            if middlerow(i)<=0.5 && middlerow(i+1)>=0.5
                startCircleInd=i;
            end

            if middlerow(i)>=0.5 && middlerow(i+1)<=0.5
                endCircleInd=i;
            end
        end
        
        disp(['The time is ' num2str(k*dt)])
        startCircleValue=(0.5-middlerow(startCircleInd))/(middlerow(startCircleInd+1)-middlerow(startCircleInd))+startCircleInd;
        endCircleValue=(0.5-middlerow(endCircleInd))/(middlerow(endCircleInd+1)-middlerow(endCircleInd))+endCircleInd;

        disp(['phi_1 circle size is ' num2str(endCircleValue-startCircleValue)])
        aa=[aa; k*dt endCircleValue-startCircleValue];
    
    end
end






