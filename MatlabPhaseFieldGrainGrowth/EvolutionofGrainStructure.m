% Evolution of grain structures

clc; clear; close all;

for gg=[2:4 8 10]
GRAINS=gg;

subfldr={['/GG_Q3/' num2str(GRAINS) 'GRAINS/']};
sf_dr=subfldr{1};

dt=0.01;
gs=1; %grid size
timesteps=5e4;
plotMOD=5e2;
MG=200; %number of cells to divide width/height

middle_grain=floor(MG/2);

[xm, ym]=meshgrid(-MG/2:gs:MG/2); %phi locations

for g=1:GRAINS
phi(:,:,g)=(rand(MG+1,MG+1)-0.5)/1500;
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
% ff1=figure('Visible','off');
% contourf(xm,ym,flipud(reshape(philong(:,1)', MG+1, MG+1)'), 'linestyle','none'); colorbar;
% caxis([0 1])
% title(['Time is: ' num2str(0)],'Interpreter','latex','FontSize',16)
% xlabel('X','Interpreter','latex','FontSize',16)
% ylabel('Y','Interpreter','latex','FontSize',16)
% axis equal;
% saveas(ff1,[pwd sf_dr sprintf('FIG%d_G1.png',scount)]);
% savefig(ff1,[pwd sf_dr sprintf('FIG%d_G1.fig',scount)]);
% ff2=figure('Visible','off');
% contourf(xm,ym,flipud(reshape(philong(:,2)', MG+1, MG+1)'), 'linestyle','none'); colorbar;
% caxis([0 1])
%title(['$N=$' num2str(GRAINS) ' t: ' num2str(k*dt)],'Interpreter','latex','FontSize',16)
% xlabel('X','Interpreter','latex','FontSize',16)
% ylabel('Y','Interpreter','latex','FontSize',16)
% axis equal;
% saveas(ff2,[pwd sf_dr sprintf('FIG%d_G2.png',scount)]);
% savefig(ff2,[pwd sf_dr sprintf('FIG%d_G2.fig',scount)]);
scount=scount+1;

othergrainindex=zeros(GRAINS-1,GRAINS);


for i=1:GRAINS
crcTemp=circshift([1:GRAINS],GRAINS-i);
othergrainindex(:,i)=crcTemp(1:GRAINS-1)';
end

for k=1:timesteps
    
    for g=1:GRAINS
        dfdphi=-philong(:,g)+philong(:,g).^3;
        for j=1:size(othergrainindex,1)
             dfdphi=dfdphi+3*philong(:,g).*philong(:,othergrainindex(j,g)).^2;
        end
        philong(:,g)=philong(:,g)-(dt).*(dfdphi-(2).*(LapOp*philong(:,g)));
    end
    
     if mod(k,plotMOD)==0
        
        phiplot=zeros(MG+1, MG+1);
        
        for g=1:GRAINS
            phitemp=flipud(reshape(philong(:,g)', MG+1, MG+1))';
            phiplot=phiplot+phitemp.^2.*(phitemp-1).^2;
            ff=figure('Visible','off');
            contourf(xm,ym,flipud(reshape(philong(:,g)', MG+1, MG+1))', 'linestyle','none'); colorbar;
            caxis([0 1])
            title(['$\phi_' num2str(g) '$ \hspace{0.1in} t:' num2str(k*dt)],'Interpreter','latex','FontSize',16)
            xlabel('X','Interpreter','latex','FontSize',16)
            ylabel('Y','Interpreter','latex','FontSize',16)
            axis equal;
            saveas(ff,[pwd sf_dr sprintf('/phi%d/phiC_%d.png',g,scount)]);
            close all;
            
        end     
        ff=figure('Visible','off');
        contourf(xm,ym,phiplot, 'linestyle','none'); colorbar;
        caxis([0 1])
        title(['N:' num2str(GRAINS) '\hspace{0.1in} t:' num2str(k*dt)],'Interpreter','latex','FontSize',16)
        xlabel('X','Interpreter','latex','FontSize',16)
        ylabel('Y','Interpreter','latex','FontSize',16)
        axis equal;
        saveas(ff,[pwd sf_dr sprintf('GG_%d.png',scount)]);
        close all;
        
        save([ pwd sf_dr '/PhiStore' num2str(scount) '.mat'], 'philong','-v7.3');
        scount=scount+1;    
    end
end


clearvars -except gg
end


