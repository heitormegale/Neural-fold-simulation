%% example
%quieck example on how to run the simulation and what each parameter
%encodes
close all
n=0;
radius_l=[];
radius_r=[];

nneural=[40]; %number of neural cells       
mu=[7.5]; %interface energy value
m=[7.5]; %apical contractility value
n=n+1;
L_0_b=[1 1 1 1];
K=[1 1 1]*1; %spring constant vector (apical, basal, side)
g=20; %area modulus
alpha=1;
name=sprintf('m%g_b%g_k%g_mu%g_Lg%d_Li1_neural%g_reflec_scan_radius',m,g,K(3),mu,L_0_b(1),nneural);
name=erase(name,".");

NTimes=10000 %simulation time steps
Plots=1000   % when to record data e.g. every 1000 simulation time steps

[points,cells,struct_g,cells_on_i,interface_length,time1_complete,H_complete] = mechanical_simulation(K,gamma,g,m,alpha,mu,L_0_b,40+nneural,nneural,NTimes,Plots);

%
% paramets when to begin plot, how many points in time to plot and what
% size the centre of colony is

Ncell=40+nneural;
Init        = 1; %begin
NTimes      = length(points);%
StepSize    = 2; % Step between time points. 1 is minimum.
 
    
    
  
    n=0;
    H=[];
    clear F;
    
    % plotting movie, energy and radius over time
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

writerObj=VideoWriter(strcat(name,'.avi'));
writerObj.FrameRate=7;



open(writerObj);
for i=1:length(F)
    writeVideo(writerObj,F(i))
end
close(writerObj);

figure(200)
plot(1:length(H),H,'k.')
xlabel('Time');
ylabel('Energy');

figure(300)
plot(1:length(radius_l),radius_l,'k.')
xlabel('Time');
ylabel('Radius');
print(strcat(name,'_radius'),'-dtiff')
Mega_matrix.(name).Cells_on_interface=cells_on_i;
Mega_matrix.(name).Variation_cell_area=points(end).dArea;
Mega_matrix.(name).Energy_over_time=H;
Mega_matrix.(name).Length_of_interface=interface_length;
Mega_matrix.(name).r_a=points(end).q_a;
Mega_matrix.(name).r_b=points(end).q_b;
Mega_matrix.(name).points=points;
