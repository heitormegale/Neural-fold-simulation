
%---------------------------------------------------------------------
% license
%---------------------------------------------------------------------
%
% Copyright 2021 Heitor C. Megale and Eyal Kazrbrun
%
% This file is a free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This file is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% A copy of the GNU General Public License is here: <http://www.gnu.org/licenses/>.
%
%---------------------------------------------------------------------

% This code demonstartes how to run the nueral tube simulaiton function,
% and an example for choice of paramters. Running this code would produce
% an *.avi file showing time steps of the simulation and a plot of
% system energy as a funciton of simualtion time.


%%  Running the simualtion
radius_l=[];
radius_r=[];

nneural=40; % Number of neural cells
nectoderm =40;% Number of non-neural cells
Ncell=nectoderm + nneural; %ottal number of cells
mu=7.5; %Interface energy value
m=7.5; %Apical contractility value
L_0_b=[1 1 1 1]; %Rest length of each spring
K=[1 1 1]*1; %Spring constant vector (apical, basal, side)
g=20; %area modulus (for volume conservation)
name=sprintf('m%g_b%g_k%g_mu%g_Lg%d_Li1_neural%g_reflec_scan_radius',m,g,K(3),mu,L_0_b(1),nneural);
name=erase(name,".");

NTimes=10000; % Number of time steps in the simulation
Plots=1000;   % when to record data e.g. every 1000 simulation time steps

[points,cells,struct_g,cells_on_i,interface_length,time1_complete,H_complete] = mechanical_simulation(K,1,g,m,1,mu,L_0_b,Ncell,nneural,NTimes,Plots);

disp('Simulation Done');
%% Saving the movie, and plotting energy and radius over time
close all;
% paramets when to begin plot, how many points in time to plot and what
% size the centre of colony is
Init        = 1; %begin
NTimes      = length(points);%
StepSize    = 2; % Step between time points. 1 is minimum.

n=0;
H=[];
clear F;
%plotting time-points from the simulation
for t = Init : StepSize: NTimes
    n=n+1;
    v_a=points(t).q_a;
    v_b=points(t).q_b;
    H(n)=points(t).H;
    radius_r(n)=points(t).radius;
    radius_l(n)=points(t).radius;
    x0r=points(t).x_0_r;
    x0l=points(t).x_0_l;
    z0r=points(t).z_0_r;
    z0l=points(t).z_0_l;
    attachment=points(t).R;
    center=zeros(1,length(attachment));
    center(1,1:end)=(v_a(1,1:end-1)+v_a(1,2:end)+v_b(1,1:end-1)+v_b(1,2:end))/4;
    center(2,1:end)=(v_a(2,1:end-1)+v_a(2,2:end)+v_b(2,1:end-1)+v_b(2,2:end))/4;
    if t>1
        Volume(n)=0;%max(max(points(t).dArea));
    end
    center_b=zeros(1,length(attachment));
    center_b(1:end)=(v_b(1,1:end-1)+v_b(1,2:end))/2;
    ID=points(t).ID; %ID=1 -> ecad, ID=2 -> ncad
    gcf= figure(1);
    ENcell=Ncell-nneural;
    %lower
    hold on
    plot(v_b(1,1:ENcell/2+1), v_b(2,1:ENcell/2+1), '-','Color',[0.5 0.8 0.9])
    hold on   
    plot(v_b(1,ENcell/2+nneural+1:end), v_b(2,ENcell/2+nneural+1:end), '-','Color',[0.5 0.8 0.9])    
    plot(v_b(1,ENcell/2+1:ENcell/2+nneural+1), v_b(2,ENcell/2+1:ENcell/2+nneural+1), '-','Color',[0.9 0.2 0.2])
    %upper   
    hold on
    hold on
    plot(v_a(1,1:ENcell/2+1), v_a(2,1:ENcell/2+1), '-','Color',[0.5 0.8 0.9])
    hold on
    plot(v_a(1,ENcell/2+1:ENcell/2+nneural+1), v_a(2,ENcell/2+1:ENcell/2+nneural+1), '-','Color',[0.9 0.2 0.2])    
    plot(v_a(1,ENcell/2+nneural+1:end), v_a(2,ENcell/2+nneural+1:end), '-','Color',[0.5 0.8 0.9])    
    %sides
    for ii=1:ENcell/2+2
        plot([v_b(1,ii), v_a(1,ii)], [v_b(2,ii), v_a(2,ii)],'-','Color',[0.5 0.8 0.9])        
    end  
    for ii=ENcell/2+nneural+1:length(v_b)
        plot([v_b(1,ii), v_a(1,ii)], [v_b(2,ii), v_a(2,ii)],'-','Color',[0.5 0.8 0.9])       
    end
    for ii=ENcell/2+1:ENcell/2+nneural+1
        plot([v_b(1,ii), v_a(1,ii)], [v_b(2,ii), v_a(2,ii)],'-','Color',[0.9 0.2 0.2])
        
    end    
    xlim([Ncell-nneural, Ncell+nneural])
    ylim([-0.2, 50])
    titlet=sprintf('T=%g',t);
    title(titlet)
    set(gca,'xtick',[]) 
    set(gca,'ytick',[])
    box on
    hold off
    axis equal
    F(n)=getframe(gcf);
    drawnow
    clf(gcf)
end


%saving into avi file
writerObj=VideoWriter(strcat(name,'.avi'));
writerObj.FrameRate=7;

open(writerObj);
for i=1:length(F)
    writeVideo(writerObj,F(i))
end
close(writerObj);

%plotting energy as a function of time
figure
plot(1:length(H),H,'k.')
xlabel('Time');
ylabel('Energy');
saveas(gcf,'Energy_plot.pdf');

%plotting energy as a function of time
figure
plot(1:length(radius_l),radius_l,'k.')
xlabel('Time');
ylabel('Radius');
saveas(gcf,'Radius_plot.pdf');

