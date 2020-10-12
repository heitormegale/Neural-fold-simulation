%% plot squares
attachment=[];
v=[1:10]
center=[1.5:9.5]
figure(1)
attachment=center+randn(1,length(center))*0.5;
%lower
plot(v, zeros(length(v)), 'g.', 'MarkerSize',15)
hold on
plot(v, zeros(length(v)), 'k-')
%upper
plot(v, ones(length(v)), 'g.', 'MarkerSize',15)
hold on
plot(v, ones(length(v)), 'k-')
%center
plot(center, zeros(length(center)), 'r.', 'MarkerSize',15)
plot(center, 0.5*ones(length(center)), 'k.', 'MarkerSize',20)
%attachment
plot(attachment, -0.5*ones(length(center)), 'b.', 'MarkerSize',20)
%sides
for ii=1:length(v)
    plot([v(ii), v(ii)], [0, 1],'k-')
    
end
%attachement connection
for ii=1:length(center)
    plot([center(ii), attachment(ii)], [0, -0.5],'k-')
    
end
xlim([0, 11])
ylim([-1,3])
axis equal
%% real data
[points,cells] = PlasticGrowthH_newHamiltonian;
%
% paramets when to begin plot, how many points in time to plot and what
% size the centre of colony is
    Init        = 2; %begin
    NTimes      = length(points);%how many
    WindowSize  = 50;
    lmax        = 1;
    StepSize    = 2; % Step between time points. 1 is minimum.
    cmp = jet(NTimes);
    close all;
    m = zeros(1,NTimes); % mean of cell sizes
    s = zeros(1,NTimes); % std of cell sizes
for t = Init : StepSize: NTimes
    v=points(t).q;
    attachment=points(t).R;
    center=zeros(1,length(attachment));
    center(1:end)=(v(1:end-1)+v(2:end))/2;
    ID=points(t).ID; %ID=1 -> ecad, ID=2 -> ncad
  gcf= figure(1)
%lower
plot(v, zeros(length(v)), 'g.', 'MarkerSize',15)
hold on
plot(v, zeros(length(v)), 'k-')
%upper
plot(v, ones(length(v)), 'g.', 'MarkerSize',15)
hold on
plot(v, ones(length(v)), 'k-')
%center
plot(center, zeros(length(center)), 'k.', 'MarkerSize',15)
plot(center(ID==1), 0.5*ones(length(center(ID==1))), 'c.', 'MarkerSize',20)
plot(center(ID==2), 0.5*ones(length(center(ID==2))), 'r.', 'MarkerSize',20)

%attachment
plot(attachment, -0.5*ones(length(center)), 'b.', 'MarkerSize',20)
%sides
for ii=1:length(v)
    plot([v(ii), v(ii)], [0, 1],'k-')
    
end
%attachement connection
for ii=1:length(center)
    plot([center(ii), attachment(ii)], [0, -0.5],'k-')
    
end

xlim([min(points(NTimes).q)-1, max(points(NTimes).q)+1])
ylim([-1,1.5])


name=sprintf('Timepoint_%d.png',t)
saveas(gcf,name)
clf(gcf)
end


%% real data changed 2d
n=0;
    for mu=[1 2 5 10 15]
        for m=[1 2 5 10 15]
 %m=2
           %mu=5
         n=n+1;
         length_i=0.75;
    close all
    %clearvars -except K_n alpha n Volume mu m
    
% L_0_b(1)=1; %natural length of basal side epidermis glass
% L_0_b(2)=length_i; %natural length of basal side epidermis interface
% L_0_b(3)=1; %natural length of basal side neuronal interface
% L_0_b(4)=1; %natural length of basal side neuronal glass
L_0_b=[1 length_i 1 1];
K=[1 1 1]*1;
gamma=0.01;
g=5;
%m=2000;
alpha=1;
name=sprintf('m%d_b%d_k%d_mu%d_Lg%d_Li3_4',m,g,K(1),mu,L_0_b(1));
[points,cells,struct_g,cells_on_i,interface_length] = PlasticGrowthH2D_line_jump(K,gamma,g,m,alpha,mu,L_0_b);

%
% paramets when to begin plot, how many points in time to plot and what
% size the centre of colony is
    Init        = 1; %begin
    NTimes      = length(points);%how many
    WindowSize  = 50;
    lmax        = 1;
    StepSize    = 2; % Step between time points. 1 is minimum.
    cmp = jet(NTimes);
    close all;
    %m = zeros(1,NTimes); % mean of cell sizes
    s = zeros(1,NTimes); % std of cell sizes
    n=0;
    clear F;
    H=[];
for t = Init : StepSize: NTimes
    n=n+1;
    v_a=points(t).q_a;
    v_b=points(t).q_b;
    H(n)=points(t).H;
    radius=points(t).radius;
    x0r=points(t).x_0_r;
    x0l=points(t).x_0_l;
    %rmsA(n)=points(t).rmsA;
    attachment=points(t).R;
    center=zeros(1,length(attachment));
    center(1,1:end)=(v_a(1,1:end-1)+v_a(1,2:end)+v_b(1,1:end-1)+v_b(1,2:end))/4;
    center(2,1:end)=(v_a(2,1:end-1)+v_a(2,2:end)+v_b(2,1:end-1)+v_b(2,2:end))/4;
    
    center_b=zeros(1,length(attachment));
    center_b(1:end)=(v_b(1,1:end-1)+v_b(1,2:end))/2;
    ID=points(t).ID; %ID=1 -> ecad, ID=2 -> ncad
  gcf= figure(1)
%lower
plot(v_b(1,:), v_b(2,:), 'g.', 'MarkerSize',15)
hold on
plot(v_b(1,:), v_b(2,:), 'k-')
%upper
plot(v_a(1,:), v_a(2,:), 'g.', 'MarkerSize',15)
hold on
plot(v_a(1,:), v_a(2,:), 'k-')
%center
%plot(center_b, zeros(length(center)), 'k.', 'MarkerSize',15)
plot(center(1,ID==1), center(2,ID==1), 'c.', 'MarkerSize',20)
plot(center(1,ID==2), center(2,ID==2), 'r.', 'MarkerSize',20)

%attachment
%plot(attachment, -0.1*ones(length(center)), 'b.', 'MarkerSize',20)
%sides
for ii=1:length(v_a)
    plot([v_b(1,ii), v_a(1,ii)], [v_b(2,ii), v_a(2,ii)],'k-')
    
end
%attachement connection
for ii=1:length(center)
    %plot([center_b(ii), attachment(ii)], [0, -0.1],'k-')
    
end
%circle
circle(x0r,0,radius)   
 
circle2(x0l,0,radius) 
xlim([9, 51])
ylim([-0.2,10])
hold off
F(n)=getframe(gcf);
drawnow
% name=sprintf('Timepoint_%d.png',t)
% saveas(gcf,name)
% clf(gcf)
end

writerObj=VideoWriter(strcat(name,'.avi'));
writerObj.FrameRate=7;

%Volume(n)=max(max(points(401).dArea));

open(writerObj);
for i=1:length(F)
    writeVideo(writerObj,F(i))
end
close(writerObj);

figure(200)
plot(1:length(H),H,'k.')
Mega_matrix.(name).Cells_on_interface=cells_on_i;
Mega_matrix.(name).Variation_cell_area=points(401).dArea;
Mega_matrix.(name).Energy_over_time=H;
Mega_matrix.(name).Length_of_interface=interface_length;
      end
  end

%figure(400)
%plot(1:length(rmsA),rmsA,'k.')
%%
figure(500)
n=0;
y=1.0;
     for mu=[1 2 5 10 15]
        for m=[1 2 5 10 15]
 length_i=1;        
         % L_0_b(1)=1; %natural length of basal side epidermis glass
% L_0_b(2)=length_i; %natural length of basal side epidermis interface
% L_0_b(3)=1; %natural length of basal side neuronal interface
% L_0_b(4)=1; %natural length of basal side neuronal glass
L_0_b=[1 length_i 1 1];
K=[1 1 1]*1;
gamma=0.01;
g=5;
        
       

        name=sprintf('m%d_b%d_k%d_mu%d_Lg%d_Li3_4.avi',m,g,K(1),mu,L_0_b(1));
        n=n+1;
        h=subplot(5,5,n);
       %set(h,'XTick',[],'YTick',[],'TickLength',[0,0], 'Position', [x y .1 .15]);
        
        obj = VideoReader(name);
        video = obj.read();
        
        imshow(video(:,:,:,200))
        title(sprintf('Miosin=%d jump=%d',m,mu))
        hold on
     end
     end
         
saveas(figure(500),'table.tif')
%% phase space
color=[0 0 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0; 0 0 0;0 1 0;0 1 0;1 0 0;1 0 0; 0 0 0; 0 1 0; 0 1 0; 0 1 0; 1 0 0; 0 0 0; 0 1 0; 1 0 0; 1 0 0; 1 0 0; 0 0 0; 0 1 0; 1 0 0; 1 0 0; 1 0 0]
color = [ color; 0 1 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0]
color=[ color; 0 1 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0;0 1 0; 0 1 0; 0 1 0; 0 1 0; 1 0 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0; 1 0 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0; 1 0 0]
n=0
for ll_i=[1 0.5 0.75]
for mu=[1 2 5 10 15]
        for m=[1 2 5 10 15]
       n=n+1     
figure(22)
plot3(m,mu,ll_i,'.','MarkerFaceColor',color(n,:),'MarkerEdgeColor',color(n,:))
hold on
end
end
end

%% real data changed 
[points,cells] = PlasticGrowthH_newHamiltonian;
%
% paramets when to begin plot, how many points in time to plot and what
% size the centre of colony is
    Init        = 1; %begin
    NTimes      = length(points);%how many
    WindowSize  = 50;
    lmax        = 1;
    StepSize    = 2; % Step between time points. 1 is minimum.
    cmp = jet(NTimes);
    close all;
    m = zeros(1,NTimes); % mean of cell sizes
    s = zeros(1,NTimes); % std of cell sizes
for t = Init : StepSize: NTimes
    v_a=points(t).q_a;
    v_b=points(t).q_b;
    attachment=points(t).R;
    center=zeros(1,length(attachment));
    center(1:end)=(v_a(1:end-1)+v_a(2:end)+v_b(1:end-1)+v_b(2:end))/4;
    
    center_b=zeros(1,length(attachment));
    center_b(1:end)=(v_b(1:end-1)+v_b(2:end))/2;
    ID=points(t).ID; %ID=1 -> ecad, ID=2 -> ncad
  gcf= figure(1)
%lower
plot(v_b(:), zeros(length(v_b)), 'g.', 'MarkerSize',15)
hold on
plot(v_b(:), zeros(length(v_b)), 'k-')
%upper
plot(v_a(:), ones(length(v_a)), 'g.', 'MarkerSize',15)
hold on
plot(v_a(:), ones(length(v_a)), 'k-')
%center
%plot(center_b, zeros(length(center)), 'k.', 'MarkerSize',15)
plot(center(ID==1), 0.5*ones(length(center(ID==1))), 'c.', 'MarkerSize',20)
plot(center(ID==2), 0.5*ones(length(center(ID==2))), 'r.', 'MarkerSize',20)

%attachment
plot(attachment, -0.1*ones(length(center)), 'b.', 'MarkerSize',20)
%sides
for ii=1:length(v_a)
    plot([v_b(ii), v_a(ii)], [0, 1],'k-')
    
end
%attachement connection
for ii=1:length(center)
    plot([center_b(ii), attachment(ii)], [0, -0.1],'k-')
    
end

%xlim([min(points(NTimes).q_b(1,:))-2, max(points(NTimes).q_b(1,:))+2])
%ylim([-0.2,1.7])


name=sprintf('Timepoint_%d.png',t)
saveas(gcf,name)
clf(gcf)
end
%%
kappa=0.5;
r=1:100
    N = length(r);

    % inner part
    i = 2:N-1;
    j = 2:N-1;
    s = ones(1,N-2)*( 2 + kappa/2); 
    
    i = [i (2:N-1)-1];
    j = [j (2:N-1)];
    s = [s ones(1,N-2)*(-1 + kappa/4)];

    i = [i (2:N-1)+1];
    j = [j (2:N-1)];
    s = [s ones(1,N-2)*(-1 + kappa/4)];
    
    % boundary
    i = [i 1];
    j = [j 1];
    s = [s  1 + kappa/4];

    i = [i 2];
    j = [j 1];
    s = [s -1 + kappa/4];

    i = [i N];
    j = [j N];
    s = [s  1 + kappa/4];

    i = [i N-1];
    j = [j N];
    s = [s -1 + kappa/4];    
    dhs = sparse(i,j,s,N,N);    
   %%
   dhl = zeros(length(ll)+1,1);    
    dhl(1:end-1) = kappa*0.5*rr;
    dhl(2:end)   = dhl(2:end) + kappa*0.5*rr;
    dhl(1:end-1) = dhl(1:end-1) - ll;
    dhl(2:end)   = dhl(2:end)   + ll;
    
    %external forces
    dhl(1)   = dhl(1)   - fext;
    dhl(end) = dhl(end) + fext;
    
function h = circle(x,y,r)
hold on
th = 0:pi/50:pi/2;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold on
end
function h = circle2(x,y,r)
hold on
th = pi/2:pi/50:pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold on
end
