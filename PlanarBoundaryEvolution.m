% Simulates the phase field equations to show Planar Boundary Evolution between Two Grains

clc; clear; close all;

subfldr={'/GG_Q1/'};
sf_dr=subfldr{1};

dt=0.01;
timesteps=1e6;
plotMOD=5e4;
MG=200; %number of cells to divide width/height

GRAINS=2;

gs=1; %grid size

[xm, ym]=meshgrid(-MG/2:gs:MG/2); %phi locations

phi(:,:,1)=double(xm>0);
phi(:,:,2)=double(xm<=0);

philong=zeros((MG+1)^2, GRAINS);

for g=1:GRAINS
    philong(:,g)=flipud(reshape(flipud(phi(:,:,g)'),(MG+1)^2,1));
end

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
ff=figure('Visible','off');
plot(xm(MG+1,:),philong(1:MG+1,1),'.-r', 'LineWidth', 1)
hold on
plot(xm(MG+1,:),philong(1:MG+1,2),'.-b','LineWidth', 1)
leg1 = legend('$\phi_1$','$\phi_2$','Location','southeast');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17);
title(['Time is: ' num2str(0)],'Interpreter','latex','FontSize',16)
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('$\phi$','Interpreter','latex','FontSize',16)
axis([-MG/2 MG/2 -0.2 1.2])
saveas(ff,[pwd sf_dr sprintf('FIG%d.png',scount)]);
savefig(ff,[pwd sf_dr sprintf('FIG%d.fig',scount)]);
scount=scount+1;

othergrainindex=[2 1];

for k=1:timesteps
    
    for g=1:GRAINS
        dfdphi=-philong(:,g)+philong(:,g).^3+3*philong(:,g).*philong(:,othergrainindex(g)).^2;
        philong(:,g)=philong(:,g)-(dt).*(dfdphi-(2).*(LapOp*philong(:,g)));
    end
    
    if mod(k,plotMOD)==0
        ff=figure('Visible','off');
        plot(xm(MG+1,:),philong(1:MG+1,1),'.-r', 'LineWidth', 1)
        hold on
        plot(xm(MG+1,:),philong(1:MG+1,2),'.-b','LineWidth', 1)
        leg1 = legend('$\phi_1$','$\phi_2$','Location','southeast');
        set(leg1,'Interpreter','latex');
        set(leg1,'FontSize',17);
        title(['Time is: ' num2str(k*dt)],'Interpreter','latex','FontSize',16)
        xlabel('X','Interpreter','latex','FontSize',16)
        ylabel('$\phi$','Interpreter','latex','FontSize',16)
        axis([-MG/2 MG/2 -0.2 1.2])
        saveas(ff,[pwd sf_dr sprintf('FIG%d.png',scount)]);
        savefig(ff,[pwd sf_dr sprintf('FIG%d.fig',scount)]);
        save([ pwd sf_dr '/PhiStore' num2str(scount) '.mat'], 'philong','-v7.3');
        scount=scount+1;
        
    end
end






