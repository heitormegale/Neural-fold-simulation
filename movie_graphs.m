%% generate plots


%% real data changed 2d z_0
%initial=points_original(10);
n=0;
radius_l=[];
radius_r=[];
 for nneural=[40]        
  mu=[7.5];%,0.5,7.5,0.1];%[15 20]%[0.01 0.1 0.5 1 2.5]
        m=[7.5];%,%0.25,1.25,18]%,7.5,7.5,0.1];%[5 7.5 10 15 20]


n=n+1;
length_i=1;
close all
L_0_b=[1 length_i 1 1];
K=[1 1 1]*1;
gamma=0.01;
g=20;
NTimes=4800000;
Plots=1000;
%m=2000;
alpha=1;
nneural=nneural;
name=sprintf('m%g_b%g_k%g_mu%g_Lg%d_Li1_neural%g_reflec_scan_radius',m,g,K(3),mu,L_0_b(1),nneural);
name=erase(name,".");
[points,cells,struct_g,cells_on_i,interface_length,time1_complete,H_complete] = mechanical_simulation(K,gamma,g,m,alpha,mu,L_0_b,40+nneural,nneural,NTimes,Plots);

%
% paramets when to begin plot, how many points in time to plot and what
% size the centre of colony is

Ncell=40+nneural;
    Init        = 1; %begin
    NTimes      = length(points);%how manyradius*sin(L(jj)/radius) +z_0_l
    WindowSize  = 50;
    lmax        = 1;
    StepSize    = 2; % Step between time points. 1 is minimum.
    cmp = jet(NTimes);
    close all;
    
    s = zeros(1,NTimes); % std of cell sizes
    n=0;
    clear F;
    H=[];
    Volume=[];
    Volume(1)=0;
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
Mega_matrix.(name).Variation_cell_area=points(401).dArea;
Mega_matrix.(name).Energy_over_time=H;
Mega_matrix.(name).Length_of_interface=interface_length;
Mega_matrix.(name).r_a=points(end).q_a;
Mega_matrix.(name).r_b=points(end).q_b;
Mega_matrix.(name).points=points;

end
 
%% energy
name=sprintf('m%g_b%g_k%g_mu%g_Lg%d_Li1_neural%g_reflec_scan',7.5,20,1,7.5,1,30);
            name=erase(name,".");
            for nn=1:1601
               energy(nn)= Super_Mega_matrix.(name).points(nn).H;
               time(nn)=nn;
               radius(nn)=Super_Mega_matrix.(name).points(nn).radius;
               if nn~=1
               theta(:,nn)=Super_Mega_matrix.(name).points(nn).theta(1:36);
               n_cell(nn)=length(find((~isnan(theta(:,nn)))==1)); 
               else
                   
                   n_cell(nn)=0;
               end
            end
            figure(44)
            plot(time,energy,'k.','MarkerSize',10)
            ax=gca;
            ax.FontSize = 14
            xlabel('Time','FontSize', 22)
            ylabel('Energy','FontSize', 22)
            
            box on
            grid on
            xlim([0, 1200])
            figure(45)
            plot(time,radius,'k.','MarkerSize',10)
            ax=gca;
            ax.FontSize = 14
            xlabel('Time','FontSize', 22)
            ylabel('Interface Radius','FontSize', 22)
            
            box on
            grid on
            ylim([0, 16])
            xlim([0, 1200])
            figure(46)
            plot(time,n_cell,'k.','MarkerSize',10)
            ax=gca;
            ax.FontSize = 14
            xlabel('Time','FontSize', 22)
            ylabel({'Interface Length'; '(number of cells)'},'FontSize', 22)
           
            box on
            grid on
            xlim([0, 1200])
            ylim([0, 20])
            z=[];

%% curvature and others
    m_vec=[7.5, 0.1,7.5];
    mu_vec=[7.5, 2.5 ,0.5];
    time=[1600, 1600, 1600]
    for n=1:length(m_vec)
        m=m_vec(n);
        mu=mu_vec(n);
        L_0_b=[1 1 1 1];
        nneural=30;

        K=[1 1 1]*1;
        gamma=0.01;
        b=20;
      
        name=sprintf('m%g_b%g_k%g_mu%g_Lg%d_Li1_neural%g_reflec_scan',m,b,K(3),mu,L_0_b(1),nneural);
        name=erase(name,".");
        apical=Super_Mega_matrix.(name).points(time(n)).q_a;
        basal=Super_Mega_matrix.(name).points(time(n)).q_b;   
        theta=Super_Mega_matrix.(name).points(time(n)).theta;
        radius=Super_Mega_matrix.(name).points(time(n)).radius;
        
        % interface length
        dist_1=sqrt((apical(1,21)-basal(1,21))^2+(apical(2,21)-basal(2,21))^2);
        dist_2=(theta(21)-min(theta))*radius;
        inter_length(n)=sum([dist_1 dist_2],'omitnan');
        inter_length=[17, 13,2]
        
        % neural area
        side_length=[];
        for jj=21:50
            side_length=[ side_length, sqrt((apical(1,jj)-apical(1,jj+1))^2+(apical(2,jj)-apical(2,jj+1))^2)];
        end
        
        neural_area(n)=sum(side_length);
        
        
        % curvature
        K_vec=[]
        for jj=21:49
            x1=apical(1,jj);
            x2=apical(1,jj+1);
            x3=apical(1,jj+2);
            y1=apical(2,jj);
            y2=apical(2,jj+1);
            y3=apical(2,jj+2);
         K = 2*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)) ./ ...
         sqrt(((x2-x1).^2+(y2-y1).^2)*((x3-x1).^2+(y3-y1).^2)*((x3-x2).^2+(y3-y2).^2));
        K_vec=[K_vec, K];
        end
        curvature(n)=mean(K_vec);
    end
    cmap = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560];
%area
% figure(7)
set(gcf,'color','w','units','centimeters');
set(gcf,'paperposition',[1 1 9 8]);
b=bar(neural_area,'facecolor','flat');
b.CData(1,:) = cmap(1,:);
b.CData(2,:) = cmap(2,:);
b.CData(3,:) = cmap(3,:);
hold on
set(gca,'xticklabel',[])
box on;
 set(gca,'fontsize',20);
 ylabel({'Neural Area'});
 set(gca,'xtick',[1 2 3],'xticklabel',{'CTRL','Low \alpha' ,'Low \epsilon'});
 
 
%curvature
 gcf=figure(8)
set(gcf,'color','w','units','centimeters');
set(gcf,'paperposition',[1 1 9 8]);
b=bar(curvature,'facecolor','flat');
b.CData(1,:) = cmap(1,:);
b.CData(2,:) = cmap(2,:);
b.CData(3,:) = cmap(3,:);
hold on
set(gca,'xticklabel',[])
set(gca, 'YScale', 'log')
box on;
 set(gca,'fontsize',20);
 ylabel({'Curvature'});
 set(gca,'xtick',[1 2 3],'xticklabel',{'CTRL','L. Contrac.' ,'L. Interface'});
 
 %length
figure(10)
set(gcf,'color','w','units','centimeters');
set(gcf,'paperposition',[1 1 9 8]);
b=bar(inter_length,'facecolor','flat');
b.CData(1,:) = cmap(1,:);
b.CData(2,:) = cmap(2,:);
b.CData(3,:) = cmap(3,:);
hold on
set(gca,'xticklabel',[])
box on;
set(gca,'fontsize',20);
ylabel({'Interface Length'});
 %set(gca,'xtick',[1 2 3],'xticklabel',{'CTRL','L. Contrac.' ,'L. Interface'});
 
 figure(16)
set(gcf,'color','w','units','centimeters');
set(gcf,'paperposition',[10 10 20 5]);

clear h
b=bar([1 2],[ 1 2; 1 2;1 2]);

l=legend({'CTRL ';'Low Contractility';' Low Interface Energy '});
l.Orientation='horizontal';
set(gca,'fontsize',40);
axis([3 10 0 1]);
axis off
%% phase space
figure(33)
p1=loglog(7.5,7.5,'.','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerSize',10)
hold on
p2=loglog(1.25,7.5,'.','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',10)
p3=loglog(0.1,0.1,'.','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor',[0.3 0.3 0.3],'MarkerSize',10)
%p4=loglog(12.5,0.9,'.','MarkerFaceColor',[0.13 0.7 0.13],'MarkerEdgeColor',[0.13 0.7 0.13],'MarkerSize',10)
length_i=1;
hold on
for mu=[0.01 0.1 0.5 1 2.5 5 7.5 10 15 20]
        for m=[0.01 0.1 0.5 1 2.5 5 7.5 10 15 20]
            
              
            L_0_b=[1 1 1 1];
            nneural=30;

            K=[1 1 1]*1;
            gamma=0.01;
            b=20;
      
            name=sprintf('m%g_b%g_k%g_mu%g_Lg%d_Li1_neural%g_reflec_scan',m,b,K(3),mu,L_0_b(1),nneural);
            name=erase(name,".");  
            
            if isfield(Super_Mega_matrix,name)
                position=Super_Mega_matrix.(name).points(end).q_a;
                z_0=Super_Mega_matrix.(name).points(end).z_0_l;
                radius=Super_Mega_matrix.(name).points(end).radius;
                z=Super_Mega_matrix.(name).points(end).q_b(2,:);
            else  
                if isfield(Mega_matrix_2,name)
                    position=Mega_matrix_2.(name).points(end).q_a;
                    z_0=Mega_matrix_2.(name).points(end).z_0_l;
                    radius=Mega_matrix_2.(name).points(end).radius;
                else
                    
                   position=zeros(2,71);
                   z_0=0;
                   radius=0;
                   name
                end
            end
            
            if position(1,51)-position(1,21)<1
                
                color=[1 0 0];
            else
            
                color=[0.3 0.3 0.3];
            end
            if position(1,21)==0
                color=[0.3 0.3 0.3];
            end
            
            if z_0>0.9*radius
                
               color=[0.3 0.3 0.3]; 
            end
            
            if z(20)>0 && z(22)>0
                
              color=[0 0 1];  
            end
            if position(1,51)-position(1,21)<1 && position(1,21)~=0
                
                color=[1 0 0];
            end
            figure(33)
            loglog(m,mu,'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',10)
            hold on
            xlabel('Myosin')
            ylabel('Interface')
            xlim([0.05 25])
            set(gca, 'XTick', [0.01 0.1 1 2.5 5 10 20])
            ylim([0.05 25])
            set(gca, 'YTick', [0.01 0.1 1 2.5 5 10 20])
        end
end


m_vec=[0.75,0.75,1,1.3,1.3,1.3,2.5,5,7.5,8.5,10,15,20,1.75,1.75,1.75,3.5,3.5,6,7.5,12.5,10,10, 0.9,0.9,0.9,1.15,1.15,1.15,1.75,1.75, 1.75,15,12.5,15,12.5,20,15,12.5,15,20];
mu_vec=[20,15,12.5,10,7.5,5,3.5,3.5,1.3,1,0.75,0.75,0.75,10,7.5,6,5,3.5,2.5,2,1,1.3,2, 20,15,12.5,20,15,12.5,12.5,15,20,0.9,0.9,1.15,1.15,1.15,2,1.75,1.75,1.75];
          for n=1:length(m_vec)
              mu=mu_vec(n);
              m=m_vec(n);
            L_0_b=[1 length_i 1 1];

            K=[1 1 1]*1;
            gamma=0.01;
            b=20;
      
            name=sprintf('m%g_b%g_k%g_mu%g_Lg%d_Li1_neural%g_reflec_scan',m,b,K(3),mu,L_0_b(1),nneural);
            name=erase(name,".");  
            
            if isfield(Super_Mega_matrix,name)
                position=Super_Mega_matrix.(name).points(end).q_a;
                z_0=Super_Mega_matrix.(name).points(end).z_0_l;
                radius=Super_Mega_matrix.(name).points(end).radius;
                z=Super_Mega_matrix.(name).points(end).q_b(2,:);
            else  
                if isfield(Mega_matrix_4,name)
                    position=Mega_matrix_4.(name).points(end).q_a;
                    z_0=Mega_matrix_4.(name).points(end).z_0_l;
                    radius=Mega_matrix_4.(name).points(end).radius;
                else
                    
                   position=zeros(2,71); 
                   z_0=0;
                   radius=0;
                   a=1+a
                end
            end
            
            if position(1,51)-position(1,21)<1
                
                color=[1 0 0];
            else
            
                color=[0.3 0.3 0.3];
            end
            if position(1,21)==0
                color=[0.3 0.3 0.3];
            end
            if z_0>0.9*radius
                
               color=[0.3 0.3 0.3]; 
            end
            if z(20)>0 && z(22)>0
                
              color=[0 0 1];  
            end
            if position(1,51)-position(1,21)<1 && position(1,21)~=0
                
                color=[1 0 0];
            end
            figure(33)
            loglog(m,mu,'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',10)
            hold on
            xlabel('\alpha')
            ylabel('\epsilon')
            xlim([0.05 25])
            set(gca, 'XTick', [0.01 0.1 1 2.5 5 10 20])
            ylim([0.05 25])
            set(gca, 'YTick', [0.01 0.1 1 2.5 5 10 20])
          end
          
          
for m=[0.1, 0.5, 1, 1.3, 1.75, 2.5, 5, 7.5, 10, 15, 20, 12,14,17,6,3,4,0.25,1.25,18]
    mu=7.5;
    L_0_b=[1 length_i 1 1];

            K=[1 1 1]*1;
            gamma=0.01;
            b=20;
      
            name=sprintf('m%g_b%g_k%g_mu%g_Lg%d_Li1_neural%g_reflec_scan',m,b,K(3),mu,L_0_b(1),nneural);
            name=erase(name,".");  
            
            if isfield(Super_Mega_matrix,name)
                position=Super_Mega_matrix.(name).points(end).q_a;
                z_0=Super_Mega_matrix.(name).points(end).z_0_l;
                radius=Super_Mega_matrix.(name).points(end).radius;
                z=Super_Mega_matrix.(name).points(end).q_b(2,:);
            else  
                if isfield(Mega_matrix_4,name)
                    position=Mega_matrix_4.(name).points(end).q_a;
                    z_0=Mega_matrix_4.(name).points(end).z_0_l;
                    radius=Mega_matrix_4.(name).points(end).radius;
                else
                    
                   position=zeros(2,71); 
                   z_0=0;
                   radius=0;
                end
            end
            
            if position(1,51)-position(1,21)<1
                
                color=[1 0 0];
            else
            
                color=[0.3 0.3 0.3];
            end
            if position(1,21)==0
                color=[0.3 0.3 0.3];
            end
            if z_0>0.9*radius
                
               color=[0.3 0.3 0.3]; 
            end
            if z(20)>0 && z(22)>0
                
              color=[0 0 1];  
            end
            if position(1,51)-position(1,21)<1 && position(1,21)~=0
                
                color=[1 0 0];
            end
            figure(33)
            loglog(m,mu,'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',10)
            hold on
            xlabel('\alpha')
            ylabel('\epsilon')
            xlim([0.05 25])
            set(gca, 'XTick', [0.01 0.1 1 2.5 5 10 20])
            ylim([0.05 25])
            set(gca, 'YTick', [0.01 0.1 1 2.5 5 10 20])
          end
    
for mu=[0.1, 0.5, 1, 1.3, 2, 2.5, 5,  7.5, 10, 15, 6, 9, 14, 3, 12, 1.15, 1.25, 2.25] 
    m=7.5;
    L_0_b=[1 length_i 1 1];

            K=[1 1 1]*1;
            gamma=0.01;
            b=20;
      
            name=sprintf('m%g_b%g_k%g_mu%g_Lg%d_Li1_neural%g_reflec_scan',m,b,K(3),mu,L_0_b(1),nneural);
            name=erase(name,".");  
            
            if isfield(Super_Mega_matrix,name)
                position=Super_Mega_matrix.(name).points(end).q_a;
                z_0=Super_Mega_matrix.(name).points(end).z_0_l;
                radius=Super_Mega_matrix.(name).points(end).radius;
                z=Super_Mega_matrix.(name).points(end).q_b(2,:);
            else  
                if isfield(Mega_matrix_4,name)
                    position=Mega_matrix_4.(name).points(end).q_a;
                    z_0=Mega_matrix_4.(name).points(end).z_0_l;
                    radius=Mega_matrix_4.(name).points(end).radius;
                else
                    
                   position=zeros(2,71); 
                   z_0=0;
                   radius=0;
                end
            end
            
            if position(1,51)-position(1,21)<1
                
                color=[1 0 0];
            else
            
                color=[0.3 0.3 0.3];
            end
            if position(1,21)==0
                color=[0.3 0.3 0.3];
            end
            if z_0>0.9*radius
                
               color=[0.3 0.3 0.3]; 
            end
            if z(20)>0 && z(22)>0
                
              color=[0 0 1];  
            end
            if position(1,51)-position(1,21)<1 && position(1,21)~=0
                
                color=[1 0 0];
            end
            figure(33)
            loglog(m,mu,'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',10)
            hold on
            xlabel('Apical Contractility (\alpha)','FontSize', 18)
            ylabel('Interface Energy (\epsilon)','FontSize', 18)
            xlim([0.05 25])
            set(gca, 'XTick', [0.01 0.1 1 2.5 5 10 20])
            ylim([0.05 25])
            set(gca, 'YTick', [0.01 0.1 1 2.5 5 10 20])
end
          
% manual
m_vec=[0.75, 0.75, 0.9, 0.9, 0.9, 1, 1, 1.15, 1.75, 7.5, 20, 10, 7.5, 7.5, 8.5, 10, 12.5, 15, 20, 20, 15, 20, 15, 12.5, 15, 12.5,0.1,0.5,0.1,5,0.1,0.5,1,0.25,0.1,20,20,20,15,15,15,0.5,0.5,0.9,1,1,1.25,1.75,2.5,3.5,5,7.5,10,7.5,10,12.5,0.1,0.5,1,2.5,0.5,0.1,0.25,1,2.5,5,7.5,20,20,15,15];
mu_vec=[15,20,20,15,12.5,12.5,10,20,20,20,20,1.3,1.15,1,1,1,1,1,1,0.75,0.75,1.15,1.15,1.15,0.9,0.9,20,20,15,1,1.25,1.25,1.25,2,2,1,0.75,0.5,1,0.9,0.75,15,10,12.5,10,7.5,7.5,6,3.5,3.5,2.5,1.3,1.3,1.25,1,0.9,1.3,1.3,1.3,1.3,1.75,1.75,1.75,1.75,1.75,1.75,1.75,0.5,0.75,0.75,0.9];
color_vec=[1,1,1,1,2,1,2,1,1,1,4,2,3,3,3,2,1,2,2,2,2,1,1,1,2,2,4,4,4,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,2,3,2,2,2,2,1,3,3,3,3]% 1 is red 2 is blue 3 is grey 4 is white 5 is green

for nn=1:length(m_vec)
    if color_vec(nn)==1
       color=[1 0 0] 
      loglog(m_vec(nn),mu_vec(nn),'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',10)  
    end
    if color_vec(nn)==2
       color=[0 0 1] 
      loglog(m_vec(nn),mu_vec(nn),'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',10)  
    end
    if color_vec(nn)==3
       color=[0.3 0.3 0.3] 
      loglog(m_vec(nn),mu_vec(nn),'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',10)  
    end
    if color_vec(nn)==4
       color=[1 1 1] 
      loglog(m_vec(nn),mu_vec(nn),'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',10)  
    end
    if color_vec(nn)==5
       color=[0.13 0.7 0.13] 
      loglog(m_vec(nn),mu_vec(nn),'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerSize',10)  
    end
    
end
ax=gca;
ax.FontSize = 12
loglog(7.5,7.5,'d','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerSize',5)
hold on
loglog(7.5,0.5,'d','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor',[0.3 0.3 0.3],'MarkerSize',5)
loglog(0.1,2.5,'d','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',5)
legend([p1,p2,p3],'Closure','Interface','No interface','Location','southwest')
hold off
%% radius
n=0;
m_vec=[];
radius=[];
%radius
for m=[0.1, 0.5, 1, 1.3, 1.75, 2.5, 5, 7.5, 10, 15, 20, 12,14,17,6,3,4,0.25,1.25,18]
    mu=7.5;
L_0_b=[1 1 1 1];
            nneural=30;
              n=n+1;
              m_vec(n)=m;
            K=[1 1 1]*1;
            gamma=0.01;
            b=20;
      
            name=sprintf('m%g_b%g_k%g_mu%g_Lg%d_Li1_neural%g_reflec_scan',m,b,K(3),mu,L_0_b(1),nneural);
            name=erase(name,".");  
            radius(n)=15;
            if isfield(Super_Mega_matrix,name)
                
                radius(n)=Super_Mega_matrix.(name).points(end).radius;
                
            end

            figure(22)
           plot(m,radius(n),'k.','MarkerSize',10)
            hold on
            xlabel('\alpha')
            ylabel('Radius')
            
end
[xData, yData] = prepareCurveData( m_vec, radius );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x+d)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [15 1 5 -1];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

plot( fitresult,'k--');
hold on
plot(m_vec,radius,'k.','MarkerSize',20)
           
legend off
xlabel('Apical contractility (\alpha)','FontSize', 18)
ylabel('Interface Radius','FontSize', 18)
ylim([0, 16])
ax=gca;
ax.FontSize = 14
box on
grid on
%% number of cells
n=0;
mu_vec=[];
n_cell=[];
for mu=[0.1, 0.5, 1, 1.3, 2, 2.5, 5,  7.5, 10, 15, 6, 9, 14, 3, 12, 1.15, 1.25, 2.25] %[0.1, 0.5,0.75,0.9, 1,1.15,1.75,2, 2.5, 5,  7.5, 10, 15]%
    m=7.5;
    
L_0_b=[1 1 1 1];
            nneural=30;
            n=n+1;
            K=[1 1 1]*1;
            gamma=0.01;
            b=20;
            mu_vec(n)=mu;
            name=sprintf('m%g_b%g_k%g_mu%g_Lg%d_Li1_neural%g_reflec_scan',m,b,K(3),mu,L_0_b(1),nneural);
            name=erase(name,".");  
            z=0
            if isfield(Super_Mega_matrix,name)
                
                z=Super_Mega_matrix.(name).points(end).q_b(2,:);
                theta=Super_Mega_matrix.(name).points(end).theta(1:36);
                
            end

            n_cell(n)=length(find((~isnan(theta))==1));
            figure(24)
            
            
            
end
[xData, yData] = prepareCurveData( mu_vec, n_cell );
ft = fittype( 'a*exp(-b*x+d)+c', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', [1 2] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [-1 0.5 20 0.0577829519793382];
opts.Upper = [0 1.3 Inf Inf];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.

plot( fitresult,'k--');
hold on
plot(mu_vec,n_cell,'k.','MarkerSize',20)
           
legend off
xlabel('Interface Energy (\epsilon)','FontSize', 18)
ylabel('Interface Size (number of cells)','FontSize', 18)
ylim([0, 25])
ax=gca;
ax.FontSize = 14
box on
grid on
