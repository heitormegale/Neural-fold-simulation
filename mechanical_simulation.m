function [points,cells,struct_g,cells_on_i,interface_length,time1_complete,H] = mechanical_simulation(K,gamma,g,m,alpha,mu,L_0_b,Ncell,Nnumber,NTimes,Plots)

    
    
    % spring mechancis model for tissue/colony growth
    % r     positions of cells
    % rr    attachment points
    % ll     length of a cell
    % kappa attachment strength

    
    PlasticGrowth = 0; % Growth in response to stress = 1; Constant Growth = 0;
    
    % parameters:

       
        NTimes = NTimes; % Number of simulation runs;
        Plots  = Plots; % When to plot distribution
        
        Enumber=Ncell-Nnumber;
        Nnumber=Nnumber;
        Ncells = Ncell+1;   % number of edges at start
        sigma  = 1;    % standard deviation of initial condition
        kappa  = 1;    % attachment strength
        gamma0 = 0.04;   % cell division rate 
        eta    = 0;    % cell growth rate %0.08 %0 is set for no cell growth
        nu     = 1;   % friction coefficient 
        dtrr   = 0.01; % time step for attachment point relaxation
        fext   = 0;    % external force on cell 1 and N %1
        k_a    =K(1); %spring constant of apical side springs
        k_b    =K(2); %spring constant of basal side springs
        k_s    = K(3); %spring constant of side springs
        gamma  =gamma;
        lambda = 0.01;
        g = g; %volume stiffnes
        delta_x=0.005;
        lmax   = 2;     % upper threshold for division
        lmin   = .35;    % lower threshold for division
        m=m; %miosin force
        radius=15;
        x_0=(Ncells-1);
        x_0_l=x_0-Nnumber/2+radius;
        x_0_r=x_0+Nnumber/2-radius;
        z_0=0;
        z_offset_l=0;
        z_offset_r=0;
         x_offset_l=0;
        x_offset_r=0;
        alpha=alpha;
        mu=mu; %interface constant
        ll_b_e_g=L_0_b(1); %natural length of basal side epidermis glass
        ll_b_e_i=L_0_b(2); %natural length of basal side epidermis interface
        ll_b_n_i=L_0_b(3); %natural length of basal side neuronal interface
        ll_b_n_g=L_0_b(4); %natural length of basal side neuronal glass
        
        time1=0;
        time2=0;
    % Chain initiation, now we also define the line in the N/E boundary
    [r_a,r_b,rr,ll_a,ll_b,ll_s,ID,A_0,step,ID_initial,ID_edge,ID_v] = ChainModelInitial(Ncells,sigma,lmax,m,ll_b_e_g,ll_b_n_g,Nnumber); % from this we get cell edges, attachement points and cell length
    

    % initial part
 x_0_l=Ncell-Nnumber/2+radius;
x_0_r=Ncell+Nnumber/2-radius;
z_0_l=0;%radius/4;
z_0_r=0;%radius/4;
radius=15%Nnumber/2;
    
    r_a(1,(Ncells-1-Nnumber)/2+1)=x_0_l-sqrt(radius^2-(1-z_0_l)^2);
    r_a(1,(Ncells-1-Nnumber)/2+Nnumber+1)=sqrt(radius^2-(1-z_0_r)^2)+x_0_r;
    
    r_b(1,(Ncells-1-Nnumber)/2+1)=x_0_l-sqrt(radius^2-(0-z_0_l)^2);
    r_b(1,(Ncells-1-Nnumber)/2+Nnumber+1)=sqrt(radius^2-(0-z_0_r)^2)+x_0_r;
    k=[k_a k_b k_s];
    
    H=Hamiltonian(ll_a,ll_b,ll_s,[r_a,r_b],kappa,rr,Ncells,k,gamma,lambda,g,A_0,step,0,radius,x_0,x_0,0,0,alpha, mu, ID_initial, A_0);
    
    points = struct('q_a',[],'R',[]);
    points(1).q_a = r_a ;
    points(1).q_b = r_b ;
    points(1).R = rr;
    points(1).l_a = ll_a;
    points(1).l_b = ll_b;
    points(1).l_s = ll_s;
    points(1).ID=ID;
    points(1).H=H;
    points(1).dH_a=[];
    points(1).dH_b=[];
    points(1).dH_a_z=[];
    points(1).dH_b_z=[];
    points(1).dArea=[];
    points(1).rmsA=0;
    points(1).radius=radius;
    points(1).x_0_r=x_0_r;
    points(1).x_0_l=x_0_l;
    points(1).z_0_r=z_0_r;
    points(1).z_0_l=z_0_l;
    points(1).interface_length=0;
    points(1).theta=0;
    points(1).theta_0_r=0;
    points(1).theta_0_l=0;
    points(1).dH_r=0;
    points(1).dH_r_total=0;
    points(1).gradnorm=0;
    points(1).Energy_ncad=[];
    points(1).Energy_ecad=[];
    points(1).Energy_left=[];
    points(1).Energy_right=[];
    points(1).apical=[];
    points(1).basal=[];
    points(1).side=[];
    points(1).area=[];
    points(1).myosin=[];
    points(1).dH_theta=[];
    points(1).delta_theta=[];
    cells = zeros(1,NTimes);
    %divposdt = [];
    cmp = jet(NTimes/Plots+1);
    idx = 1;
        x_0_r=x_0_r;
        x_0_l=x_0_l;
        z_0_r=z_0_r;
        z_0_l=z_0_l;
   

    L=zeros(length(r_a),1);
    L_a=zeros(length(r_a),1);
    L_a(ID_edge==1)=  radius*atan2(r_a(2,find(ID_edge==1))-z_0_l,-(r_a(1,find(ID_edge==1))-x_0_l));
    L_a(ID_edge==2)= radius*atan2(r_a(2,find(ID_edge==2))-z_0_r,(r_a(1,find(ID_edge==2))-x_0_r));
    
    theta=ones(length(r_a),1).*nan;
    theta_a=ones(length(r_a),1).*nan;
    theta_a(ID_edge==1)=  atan2(r_a(2,find(ID_edge==1))-z_0_l,-(r_a(1,find(ID_edge==1))-x_0_l));
    theta_a(ID_edge==2)= atan2(r_a(2,find(ID_edge==2))-z_0_r,(r_a(1,find(ID_edge==2))-x_0_r));
    theta(ID_edge==1)=atan2(r_b(2,find(ID_edge==1))-z_0_l,-(r_b(1,find(ID_edge==1))-x_0_l));
    theta(ID_edge==2)=atan2(r_b(2,find(ID_edge==2))-z_0_r,(r_b(1,find(ID_edge==2))-x_0_r));
    
    theta_a=round(theta_a,10);
    theta=round(theta,10);
    
    L(ID_edge==1)=-radius*abs(asin(z_0_l/radius))*sign(z_0_l);
    L(ID_edge==2)=-radius*abs(asin(z_0_r/radius))*sign(z_0_r);
    ok1=1;
    ok2=1;
    lim=1000;
    for i= 1:NTimes
        
       
        n=0;
        [delta_L_l, delta_L_r, ID_min_l_e, ID_min_l_n, ID_min_r_e,ID_min_r_n]=interface_length_func(L,ID_initial);
        
        if (i-time1)==lim
           ok1=1; 
        end
        if (i-time2)==lim
           ok2=1; 
        end
ok_l=1;
ok_r=1;
for jj=1:Ncells
 if jj~=length(r_a)
            x1=[r_b(1,jj) r_b(1,jj+1) r_a(1,jj+1) r_a(1,jj)];
            y1=[r_b(2,jj) r_b(2,jj+1) r_a(2,jj+1) r_a(2,jj)];
            Area(jj)=polyarea(x1,y1);
 end
end
        for jj=1:Ncells
            
            n=n+1;

            
            dH_a_x(1,n)=dH_a_func(ll_a,ll_s,r_a,r_b,kappa,jj,1,k,g,delta_x,step,A_0,Area);
            [dH_a_z(1,n)]=dH_a_func(ll_a,ll_s,r_a,r_b,kappa,jj,2,k,g,delta_x,step,A_0,Area);
            
            dH_b_x(1,n)=dH_b_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,1,k,gamma,lambda,g,delta_x,A_0,Area);
            dH_b_z(1,n)=0;%dH_b_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,2,k,gamma,lambda,g,delta_x,A_0,Area);  
            

            dH_l(1,n)=0;
            dH_r(1,n)=0;
            dH_x_0(1,n)=0;
             if ID_initial(jj)==1
                dH_b_z(1,n)=dH_b_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,2,k,gamma,lambda,g,delta_x,A_0,Area);%+dH_extra(r_b(1,jj),r_b(2,jj),radius,x_0,delta_L_l,2,jj,ID_min_l_e,ID_min_l_n,alpha);
                dH_b_x(1,n)=dH_b_x(1,n);%+dH_extra(r_b(1,jj),r_b(2,jj),radius,x_0,delta_L_l,1,jj,ID_min_l_e,ID_min_l_n,alpha);
                
                dx=-dx_dtheta(radius,theta(jj));
                dz=dz_dtheta(radius,theta(jj));
                dH_theta(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
            
                dx=-dx_dr(r_b(1,jj),r_b(2,jj),x_0_l,radius,z_0_l,L(jj),theta(jj));
                dz=dz_dr(r_b(1,jj),r_b(2,jj),x_0_l,radius,z_0_l,L(jj),theta(jj));
                dH_r(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
                
                dx=1;
                dz=0;
                dH_x_0_l(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
                
                dx=0;
                dz=1;
                dH_z_0_l(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
                if ID_edge(jj)==1
  
                    dx=-dx_dtheta(radius,theta_a(jj));
                    dz=dz_dtheta(radius,theta_a(jj));
                    dH_theta_a(1,n)=dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                    
                    dx=-dx_dr(r_a(1,jj),r_a(2,jj),x_0_l,radius,z_0_l,L_a(jj),theta_a(jj));
                    dz=dz_dr(r_a(1,jj),r_a(2,jj),x_0_l,radius,z_0_l,L_a(jj),theta_a(jj));
                    dH_r(1,n)=dH_r(1,n)+ dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                
                    dx=1;
                    dz=0;
                    dH_x_0_l(1,n)= dH_x_0_l(1,n)+ dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                
                    dx=0;
                    dz=1;
                    dH_z_0_l(1,n)=dH_z_0_l(1,n)+ dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                end
                
                    
                
            end
            if ID_initial(jj)==2
                dH_b_z(1,n)=dH_b_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,2,k,gamma,lambda,g,delta_x,A_0,Area);%+dH_extra(r_b(1,jj),r_b(2,jj),radius,x_0,delta_L_r,2,jj,ID_min_r_e,ID_min_r_n,alpha); 
                dH_b_x(1,n)=dH_b_x(1,n);%+dH_extra(r_b(1,jj),r_b(2,jj),radius,x_0,delta_L_r,1,jj,ID_min_r_e,ID_min_r_n,alpha);
                
                dx=dx_dtheta(radius,theta(jj));
                dz=dz_dtheta(radius,theta(jj));
                dH_theta(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
            
                dx=dx_dr(r_b(1,jj),r_b(2,jj),x_0_r,radius,z_0_r,L(jj),theta(jj));
                dz=dz_dr(r_b(1,jj),r_b(2,jj),x_0_r,radius,z_0_r,L(jj),theta(jj));
                dH_r(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
                
                dx=1;
                dz=0;
                dH_x_0_r(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
                
                dx=0;
                dz=1;
                dH_z_0_r(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
                if ID_edge(jj)==2
  
                    dx=dx_dtheta(radius,theta_a(jj));
                    dz=dz_dtheta(radius,theta_a(jj));
                    dH_theta_a(1,n)=dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                    
                    dx=dx_dr(r_a(1,jj),r_a(2,jj),x_0_r,radius,z_0_r,L_a(jj),theta_a(jj));
                    dz=dz_dr(r_a(1,jj),r_a(2,jj),x_0_r,radius,z_0_r,L_a(jj),theta_a(jj));
                    dH_r(1,n)=dH_r(1,n)+ dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                
                    dx=1;
                    dz=0;
                    dH_x_0_r(1,n)= dH_x_0_r(1,n)+ dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                
                    dx=0;
                    dz=1;
                    dH_z_0_r(1,n)=dH_z_0_r(1,n)+ dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                end
            end

        end

        
            if i==529000
                a=1;
            end
        N=Ncells;
        cells(i) = Ncells;
        fullgrad=[dH_a_x(1,find(ID_edge==0)) dH_a_z(1,find(ID_edge==0)) dH_b_x(1,find(ID_initial==0)) dH_b_z(1,find(ID_initial==0)) dH_theta dH_r dH_x_0_l dH_x_0_r dH_theta_a dH_z_0_l dH_z_0_r];
        gradnorm=norm(fullgrad);
        dH_r_total=sum(dH_r)/gradnorm;
        dH_x_0_l_total=sum(dH_x_0_l)/gradnorm;
        dH_x_0_r_total=sum(dH_x_0_r)/gradnorm;
        dH_z_0_l_total=sum(dH_z_0_l)/gradnorm;
        dH_z_0_r_total=sum(dH_z_0_r)/gradnorm;
        %old
        theta_0_r=-abs(asin(z_0_r/radius))*sign(z_0_r);
        theta_0_l=-abs(asin(z_0_l/radius))*sign(z_0_l);
        theta_0_r=round(theta_0_r,10);
        theta_0_l=round(theta_0_l,10);
        %rest
        
        if length(find(ID_initial(:)>0))<0
            
           dH_r_total=0; 
            
        end
        
 
        radius=radius+dtrr*(-dH_r_total)/gradnorm;
        x_offset_r=x_offset_r+dtrr*(-dH_x_0_r_total)/gradnorm;
      x_offset_l=x_offset_l+dtrr*(-dH_x_0_l_total)/gradnorm;%%here

       x_0_r=x_0_r+dtrr*(-dH_x_0_r_total)/(gradnorm);
            x_0_l=x_0_l+dtrr*(-dH_x_0_l_total)/(gradnorm);

        z_0_r=z_0_r+dtrr*(-dH_z_0_r_total)/(gradnorm);
        z_0_l=z_0_l+dtrr*(-dH_z_0_l_total)/(gradnorm);       
        struct_g(i).dH_theta=dH_theta;
        if z_0_r>=radius
            z_0_r=radius;
            ok_r=0;
        end
        if z_0_l>=radius
            z_0_l=radius;
            ok_l=0;
        end
        radius=round(radius,10);
        z_0_l=round(z_0_l,10);
        z_0_r=round(z_0_r,10);
        x_0_l=round(x_0_l,10);
        x_0_r=round(x_0_r,10);

        %new
        theta_0_r_n=-abs(asin(z_0_r/radius))*sign(z_0_r);
        theta_0_l_n=-abs(asin(z_0_l/radius))*sign(z_0_l);
       theta_0_r_n=round(theta_0_r_n,10);
       theta_0_l_n=round(theta_0_l_n,10);
     % move edges 
     
     
           
            r_a(2,1)=r_a(2,1)+dtrr*(-dH_a_z(1,1))/gradnorm;%+dArea(:,2)');
            r_a(2,end)=r_a(2,end)+dtrr*(-dH_a_z(1,end))/gradnorm;%+dArea(:,2)');
     
            if i==529000
                a=1;
            end
        for jj=2:length(r_a)-1
            
            if ID_initial(jj)==2

                theta_virt=theta(jj)+dtrr*(-dH_theta(1,jj))/gradnorm;
                virtual(2)=radius*sin(theta_virt)+ z_0_r;
                virtual(1)=radius*cos(theta_virt)+x_0_r;

                
                if theta(jj)==theta_0_r
                                if i==529000
                a=1;
            end
                    if dH_theta(jj)/gradnorm<0 && virtual(2)>0
   
                    [do_it, virtual_real]=jump(r_a,r_b,ll_s,ll_b,k,ID_edge,jj,virtual,dtrr,mu,dH_l./gradnorm,theta,theta_a);
                    a=1;
                    if do_it==1 && ok2==1 && ok_r==1
                         dH_theta(1,jj)=-(virtual_real-theta(jj))*gradnorm/dtrr;
                         
                         time2=i;
                         ok2=0;
                        if ID_edge(jj)==2
                        dH_b_x(jj-1)=-(x_0_r+sqrt(radius^2-z_0_r^2)-r_b(1,jj-1))/dtrr*gradnorm;
                        r_b(1,jj-1)=x_0_r+sqrt(radius^2-z_0_r^2);
                        ID_edge(jj-1)=1.75;
                        ID_initial(jj-1)=2;
                        theta(jj-1)=theta_0_r_n;
                        
                        dH_b_x(jj+1)=-(x_0_r+sqrt(radius^2-z_0_r^2)-r_b(1,jj+1))/dtrr*gradnorm;
                        %r_b(1,jj+1)=x_0_l-radius;
                        ID_edge(jj+1)=2.25;
                        
                        ll_b(jj-1)=ll_b_e_i;
                        ll_b(jj)=ll_b_n_i;
                        
                        end
                        if ID_edge(jj)==1.75
                        dH_b_x(jj-1)=-(x_0_r+sqrt(radius^2-z_0_l^2)-r_b(1,jj-1))/dtrr*gradnorm;
                        r_b(1,jj-1)=x_0_r+sqrt(radius^2-z_0_r^2);
                        ID_edge(jj-1)=1.75;
                        ID_initial(jj-1)=2;
                        theta(jj-1)=theta_0_r_n;
                        
                        ll_b(jj-1)=ll_b_e_i;
                        
                        end
                        if ID_edge(jj)==2.25                    
                        dH_b_x(jj+1)=-(x_0_r+sqrt(radius^2-z_0_r^2)-r_b(1,jj+1))/dtrr*gradnorm;
                        %r_b(1,jj+1)=x_0_l-radius;
                        ID_edge(jj+1)=2.25;
                        
                        
                        ll_b(jj)=ll_b_n_i;
                        end
                     else
                        
                    dH_theta(1,jj)=0;
                    theta(jj)=theta_0_r_n;
                    end
                    else

                    theta(jj)=theta_0_r_n;
                    dH_theta(jj)=0;
                    end
                end

               
                theta(jj)=theta(jj)+dtrr*(-dH_theta(1,jj))/gradnorm;
                if theta(jj)<theta_0_r_n
                    
                    theta(jj)=theta_0_r_n;
                end
             theta=round(theta,10);
                
                r_b(1,jj)= x_0_r+radius*cos(theta(jj));
                r_b(2,jj)=radius*sin(theta(jj)) +z_0_r;
                r_b(2,(r_b(2,:)<0))=0;
                
%                 if r_b(2,jj)==r_b(2,jj-1) && r_b(1,jj)==r_b(1,jj-1)
%                     
%                   r_b(2,jj)=r_b(2,jj)+0.001;
%                   r_b(1,jj)=r_b(1,jj)+0.001;
%                 end
                if ID_edge(jj)==2
                    theta_a(jj)=theta_a(jj)+dtrr*(-dH_theta_a(1,jj))/gradnorm;
                    theta_a=round(theta_a,10);
                    r_a(1,jj)=x_0_r+radius*cos(theta_a(jj));
                    r_a(2,jj)=radius*sin(theta_a(jj))+z_0_r;   
   
                else
                    r_a(1,jj)=r_a(1,jj)+dtrr*(-dH_a_x(1,jj))/gradnorm;%+dArea(:,1)');%+dtrr*noise
                    r_a(2,jj)=r_a(2,jj)+dtrr*(-dH_a_z(1,jj))/gradnorm;
                end
                
                if r_a(1,jj)<Ncell
                    r_a(1,jj)=Ncell;
                    theta_a(jj)=theta_a(jj)-dtrr*(-dH_theta_a(1,jj))/gradnorm;
                end
                 if r_b(1,jj)<Ncell
                    r_b(1,jj)=Ncell;
                    theta(jj)=theta(jj)-dtrr*(-dH_theta(1,jj))/gradnorm;
                end
            end
           
            
            
            if ID_initial(jj)==1
%                 if L(jj)==0
%                     jj=jj
%                 a=dH_l(jj)/gradnorm
%                 end
                theta_virt=theta(jj)+dtrr*(-dH_theta(1,jj))/gradnorm;
                virtual(2)=radius*sin(theta_virt)+ z_0_l;
                 virtual(1)=x_0_l-radius*cos(theta_virt);
%                 
%                contact(2)=radius*sin(theta_0_l_n)+ z_0_l;
%                contact(1)=x_0_l-radius*cos(theta_0_l_n);
%                 
%                 dist=sqrt((virtual(1)-contact(1))^2+(virtual(2)-contact(2))^2);
                
                if theta(jj)==theta_0_l 
                                if i==529000
                a=1;
            end
                    if dH_theta(jj)/gradnorm<0 && virtual(2)>0
 
                    [do_it, virtual_real]=jump(r_a,r_b,ll_s,ll_b,k,ID_edge,jj,virtual,dtrr,mu,dH_l./gradnorm,theta,theta_a);
                    
                    if do_it==1 && ok1==1 && ok_l==1
                         dH_theta(1,jj)=-(virtual_real-theta(jj))*gradnorm/dtrr;
                         
                         time1=i;
                         ok1=0;
                        if ID_edge(jj)==1
                        dH_b_x(jj-1)=-(x_0_l-sqrt(radius^2-z_0_l^2)-r_b(1,jj-1))/dtrr*gradnorm;
                        r_b(1,jj-1)=x_0_l-sqrt(radius^2-z_0_l^2);
                        ID_edge(jj-1)=0.75;
                        ID_initial(jj-1)=1;
                        theta(jj-1)=theta_0_l_n;
                        
                        dH_b_x(jj+1)=-(x_0_l-sqrt(radius^2-z_0_l^2)-r_b(1,jj+1))/dtrr*gradnorm;
                        %r_b(1,jj+1)=x_0_l-radius;
                        ID_edge(jj+1)=1.25;
                        
                        ll_b(jj-1)=ll_b_e_i;
                        ll_b(jj)=ll_b_n_i;
                        
                        end
                        if ID_edge(jj)==0.75
                        dH_b_x(jj-1)=-(x_0_l-sqrt(radius^2-z_0_l^2)-r_b(1,jj-1))/dtrr*gradnorm;
                        r_b(1,jj-1)=x_0_l-sqrt(radius^2-z_0_l^2);
                        ID_edge(jj-1)=0.75;
                        ID_initial(jj-1)=1;
                        theta(jj-1)=theta_0_l_n;
                        
                        ll_b(jj-1)=ll_b_e_i;
                        
                        end
                        if ID_edge(jj)==1.25                    
                        dH_b_x(jj+1)=-(x_0_l-sqrt(radius^2-z_0_l^2)-r_b(1,jj+1))/dtrr*gradnorm;
                        %r_b(1,jj+1)=x_0_l-radius;
                        ID_edge(jj+1)=1.25;
                        
                        
                        ll_b(jj)=ll_b_n_i;
                        end
                     else
                        
                    dH_theta(1,jj)=0;
                    theta(jj)=theta_0_l_n;
                    end
                    
                    else
                        
                    theta(jj)=theta_0_l_n;
                    dH_theta(jj)=0;  
                    end
                    
                end
%                 if L(jj) ==L_0_l && virtual(2)<=0 
%                   
%                   dH_l(1,jj)=0;
%                     L(jj)=L_0_l_n;  
%                 end
                
                theta(jj)=theta(jj)+dtrr*(-dH_theta(1,jj))/gradnorm;
                theta=round(theta,10);
                if theta(jj)<theta_0_l_n
                    
                    theta(jj)=theta_0_l_n;
                end
             
                
                r_b(1,jj)= x_0_l-radius*cos(theta(jj));
                r_b(2,jj)=radius*sin(theta(jj)) +z_0_l;
                r_b(2,(r_b(2,:)<0))=0;
                
%                 if r_b(2,jj)==r_b(2,jj-1) && r_b(1,jj)==r_b(1,jj-1)
%                     
%                   r_b(2,jj)=r_b(2,jj)+0.001;
%                   r_b(1,jj)=r_b(1,jj)+0.001;
%                 end
                if ID_edge(jj)==1
                    theta_a(jj)=theta_a(jj)+dtrr*(-dH_theta_a(1,jj))/gradnorm;
                    theta_a=round(theta_a,10);
                    r_a(1,jj)=x_0_l-radius*cos(theta_a(jj));
                    r_a(2,jj)=radius*sin(theta_a(jj))+z_0_l;   
   
                else
                    r_a(1,jj)=r_a(1,jj)+dtrr*(-dH_a_x(1,jj))/gradnorm;%+dArea(:,1)');%+dtrr*noise
                    r_a(2,jj)=r_a(2,jj)+dtrr*(-dH_a_z(1,jj))/gradnorm;
                end
                
                if r_a(1,jj)>Ncell
                    r_a(1,jj)=Ncell;
                    theta_a(jj)=theta_a(jj)-dtrr*(-dH_theta_a(1,jj))/gradnorm;
                end
                 if r_b(1,jj)>Ncell
                    r_b(1,jj)=Ncell;
                    theta(jj)=theta(jj)-dtrr*(-dH_theta(1,jj))/gradnorm;
                end
            end
            
            if ID_initial(jj)==0
            noise=.2*randn(length(r_a),1);
            r_a(1,jj)=r_a(1,jj)+dtrr*(-dH_a_x(1,jj))/gradnorm;%+dArea(:,1)');%+dtrr*noise%
            r_a(2,jj)=r_a(2,jj)+dtrr*(-dH_a_z(1,jj))/gradnorm;%+dArea(:,2)');
        if jj==50
            d=1;
        end
            noise=.2*randn(length(r_b),1);
            r_b(1,jj)=r_b(1,jj)+dtrr*(-dH_b_x(1,jj))/gradnorm;%+dArea(:,3)');%+dtrr*noise?
            r_b(2,jj)=r_b(2,jj);%+dtrr*(-dH_b_z(1,jj))/gradnorm;%+dArea(:,4)');
            r_b(2,(r_b(2,:)<0))=0;
            if 0==abs(x_0_l-sqrt(radius^2-z_0_l^2)-r_b(1,jj))
                dH_b_x(jj)=(x_0_l-sqrt(radius^2-z_0_l^2)-r_b(1,jj))/dtrr;
                ID_initial(jj)=1;
                 r_b(1,jj)=x_0_l-sqrt(radius^2-z_0_l^2);
                  theta(jj)=theta_0_l_n;
            end
            if 0==abs(x_0_r+sqrt(radius^2-z_0_r^2)-r_b(1,jj))
                dH_b_x(jj)=(x_0_r+sqrt(radius^2-z_0_r^2)-r_b(1,jj))/dtrr;
                ID_initial(jj)=2;
                r_b(1,jj)=x_0_r+sqrt(radius^2-z_0_r^2);
                theta(jj)=theta_0_r_n;
            end
            end
            if i>1
                if jj<Ncells/2 && ID_initial(jj)>0 && theta(jj)~=theta_0_l_n
               z= line_l(1)*r_a(1,jj)+line_l(2);
               if r_a(2,jj)<z-0.01
                   r_a(2,jj)=z-0.01;
               end
                end
                
                if jj>Ncells/2 && ID_initial(jj)>0 && theta(jj)~=theta_0_r_n
               z= line_r(1)*r_a(1,jj)+line_r(2);
               if r_a(2,jj)<z-0.01
                   r_a(2,jj)=z-0.01;
               end
                end
            end
%             if theta_a(jj)<-1.4
%                 theta_a(jj)=-1.4;
%             end

        end
        r_a(1,(Ncells-1)/2+1)=Ncells-1;
        r_b(1,(Ncells-1)/2+1)=Ncells-1;
        r_a=round(r_a,10);
        r_b=round(r_b,10);
        theta=round(theta,10);
        theta_a=round(theta_a,10);
%                                if (r_a(1,22)-70+(r_a(1,50)-70))~= 0%length(find(ID_edge>0))>2
%     b=1;
% end 
                 
        for jj=2:length(r_a)-1
            
            if ID_edge(jj)==2.25 
           if theta(jj)==theta(jj+1)
               theta(jj)=theta(jj)+0.0001;
           r_b(1,jj)=radius*cos(theta(jj)) + x_0_r;
           r_b(2,jj)=radius*sin(theta(jj)) + z_0_r;
           end
            end
            
            if ID_edge(jj)==1.75 
           if theta(jj)==theta(jj-1)
               theta(jj)=theta(jj)+0.0001;
           r_b(1,jj)=radius*cos(theta(jj)) + x_0_r;
           r_b(2,jj)=radius*sin(theta(jj)) + z_0_r;
           end
            end
           
            
            
            
         if ID_edge(jj)==1.25 
           if theta(jj)==theta(jj+1)
               theta(jj)=theta(jj)+0.0001;
               r_b(1,jj)= x_0_l-radius*cos(theta(jj));
               r_b(2,jj)=radius*sin(theta(jj)) +z_0_l;
           end
         end
          if ID_edge(jj)==0.75 
           if theta(jj)==theta(jj-1)
               theta(jj)=theta(jj)+0.0001;
               r_b(1,jj)= x_0_l-radius*cos(theta(jj));
               r_b(2,jj)=radius*sin(theta(jj)) +z_0_l;
           end
          end
            
        end
        r_a=round(r_a,10);
        r_b=round(r_b,10);
%                                  if (r_a(1,22)-70+(r_a(1,50)-70))~= 0%length(find(ID_edge>0))>2
%     c=1;
% end

        [line_l, line_r]=ecad_lines(r_a,ID_initial,ID_v,Ncells);
        %adhesion relaxation
        xi      = .2*randn(length(rr),1);        %Generate random noise
        xi(1)   =  -abs(xi(1))-fext;        %active traction on boundaries
        xi(end) = abs(xi(end))+fext;
        rr      = rr + dtrr  *(kappa/nu * (0.5*(r_b(1,1:end-1)+r_b(1,2:end))-rr) ); %+ sqrt(dtrr)*xi;
        %rr      = rr + dtrr  * (kappa/nu * (r-rr) ) + sqrt(dtrr)*xi; %spring from basal point to apical point above it
        
   if mod(i,10*1000) == 0
        
    reflection_a=(-r_a(:,Ncell/2:-1:1)+Ncell*2*[ones(1,Ncell/2);zeros(1,Ncell/2)]).*[ones(1,Ncell/2);-1*ones(1,Ncell/2)];
    r_a=[r_a(:,1:Ncell/2+1),reflection_a];
    reflection_b=(-r_b(:,Ncell/2:-1:1)+Ncell*2*[ones(1,Ncell/2);zeros(1,Ncell/2)]).*[ones(1,Ncell/2);-1*ones(1,Ncell/2)];
    r_b=[r_b(:,1:Ncell/2+1),reflection_b];
    
    radius=radius;
    x_0_l=x_0_l;
    x_0_r=Ncell-(x_0_l-Ncell);
    z_0_l=z_0_l;
    z_0_r=z_0_l;
    
    theta_0_r_n=theta_0_l_n;
    theta_a=[theta_a(1:Ncell/2+1);theta_a(Ncell/2:-1:1)];
    theta=[theta(1:Ncell/2+1);theta(Ncell/2:-1:1)];
    
    ID_initial=[ID_initial(1:Ncell/2+1);2*ID_initial(Ncell/2:-1:1)];
    ID_edge=[ID_edge(1:find(ID_edge==1)-1); 1; ID_edge(find(ID_edge==1)+1:Ncell/2+1); 1.75*ID_edge(Ncell/2:-1:find(ID_edge==1)+1)./1.25; 2; 2.25*ID_edge(find(ID_edge==1)-1:-1:1)./0.75 ];
            
            
   end 

            
   
        [H(i), energy, apical, basal, side, area, myosin]=Hamiltonian(ll_a,ll_b,ll_s,[r_a,r_b],kappa,rr,Ncells,k,gamma,lambda,g,A_0,step,L,radius,x_0_r,x_0_l,delta_L_l,delta_L_r,alpha, mu, ID_initial,Area);
        

        
        if i~=1

            dArea=abs(A_0-Area)./A_0*100;
        end

        time1_complete(i)=time1;
          
        if mod(i,Plots) == 0
           
            points(i/Plots+1).q_a = r_a ;
            points(i/Plots+1).q_b = r_b ;
            points(i/Plots+1).R = rr;
            points(i/Plots+1).l_a = ll_a;
            points(i/Plots+1).l_b = ll_b;
            points(i/Plots+1).l_s = ll_s;
            points(i/Plots+1).ID=ID;
            points(i/Plots+1).H=H(i);
            points(i/Plots+1).dH_a=dH_a_x;
            points(i/Plots+1).dH_b=dH_b_x;
            points(i/Plots+1).dH_a_z=dH_a_z;
            points(i/Plots+1).dH_b_z=dH_b_z;
            %points(i/Plots+1).dArea=dArea;
            points(i/Plots+1).rmsA=[];
            points(i/Plots+1).radius=radius;
            points(i/Plots+1).x_0_r=x_0_r;
            points(i/Plots+1).x_0_l=x_0_l;
            points(i/Plots+1).interface_length=max(L_a);
            points(i/Plots+1).z_0_r=z_0_r;
            points(i/Plots+1).z_0_l=z_0_l;
            points(i/Plots+1).theta=theta;
            points(i/Plots+1).theta_0_r=theta_0_r_n;
            points(i/Plots+1).theta_0_l=theta_0_l_n;
            points(i/Plots+1).dH_r=dH_r;
            points(i/Plots+1).dH_r_total=dH_r_total;
            points(i/Plots+1).gradnorm=gradnorm;
            points(i/Plots+1).Energy_ncad=energy(ID_v==2);
            points(i/Plots+1).Energy_ecad=energy(ID_v==1);
            points(i/Plots+1).Energy_left=energy(ID_initial==1);
            points(i/Plots+1).Energy_right=energy(ID_initial==2);
            points(i/Plots+1).apical=apical(ID_v==2);
            points(i/Plots+1).basal=basal(ID_v==2);
            points(i/Plots+1).side=side(ID_v==2);
            points(i/Plots+1).area=area(ID_v==2);
            points(i/Plots+1).myosin=myosin(ID_v==2);
            points(i/Plots+1).dH_theta=dH_theta;
            points(i/Plots+1).delta_theta=theta_a(21)-theta_a(51);
        end
        %if mod(i,200) == 0
        %    divposdt = [];
        %end
        dp = divpro(r_a,lmin,lmax,dtrr,gamma0);


        if sum(dp) > 0
            %ix = find(dp==1);

            %divposdt = [divposdt r(ix)'];
            [r,rr,ll,ID] = division(r,rr,ll,ID,dp);
            dhs = dhfunc(r,kappa);
            Ncells = length(r);
        end
        
    end
    %figure(4)
    %semilogy(cells)
   cells_on_i=length(find(L~=0))+2;
   interface_length=max(L_a);
end

function [line_l, line_r]=ecad_lines(r_a,ID_initial,ID_v,NCells)

ll=find(ID_initial==0); %not in interface
ll_2=find(ID_v==1); %ecad

points_id=intersect(ll,ll_2);

points_l=points_id(find(points_id<NCells/2));
x_l=r_a(1,points_l);
z_l=r_a(2,points_l);
line_l=polyfit(x_l,z_l,1);

points_r=points_id(find(points_id>NCells/2));
x_r=r_a(1,points_r);
z_r=r_a(2,points_r);
line_r=polyfit(x_r,z_r,1);

%a=NaN



end
function [do_it, virtual_real]=jump(r_a,r_b,ll_s,ll_b,k,ID_edge,jj,virtual,dtrr,mu,dH_l,theta,theta_a)
grad=dH_l(jj);
interface=-mu*1;
if ID_edge(jj)==1 || ID_edge(jj)==2
    
     virtual(1)=abs((r_a(1,jj)+1-r_b(1,jj))/3+r_b(1,jj));
 virtual(2)=abs((r_a(2,jj)+1-r_b(2,jj))/3+r_b(2,jj));
 virtual_real=(theta_a(jj)-theta(jj))/3+theta(jj);
    
%before
vector_side=[r_a(1,jj-1)-r_b(1,jj-1), r_a(2,jj-1)-r_b(2,jj-1), 0];
vector_basal=[r_b(1,jj-1)-r_b(1,jj-2), r_b(2,jj-1)-r_b(2,jj-2), 0];
vector_edge=[r_b(1,jj)-r_b(1,jj-1), r_b(2,jj)-r_b(2,jj-1), 0];
vector_inter=[r_a(1,jj)-r_b(1,jj), r_a(2,jj)-r_b(2,jj), 0];

basal1=k(2)/2*(norm(vector_basal)-ll_b(jj-2))^2;%+k(2)/2*(sqrt((r_b(1,N)-r_b(1,N-1))^2+(r_b(2,N)-r_b(2,N-1))^2)-ll_b(N-1))^2;
side1=k(3)/2*(norm(vector_side)-ll_s(jj-1))^2;
edge1=k(2)/2*(norm(vector_edge)-ll_b(jj-1))^2;
inter1=k(2)/2*(norm(vector_inter)-ll_s(jj))^2;

vector_side=[r_a(1,jj+1)-r_b(1,jj+1), r_a(2,jj+1)-r_b(2,jj+1), 0];
vector_basal=[r_b(1,jj+1)-r_b(1,jj+2), r_b(2,jj+1)-r_b(2,jj+2), 0];
vector_edge=[r_b(1,jj)-r_b(1,jj+1), r_b(2,jj)-r_b(2,jj+1), 0];

basal2=k(2)/2*(norm(vector_basal)-ll_b(jj+1))^2;%+k(2)/2*(sqrt((r_b(1,N)-r_b(1,N-1))^2+(r_b(2,N)-r_b(2,N-1))^2)-ll_b(N-1))^2;
side2=k(3)/2*(norm(vector_side)-ll_s(jj+1))^2;
edge2=k(2)/2*(norm(vector_edge)-ll_b(jj))^2;

before=basal1+side1+edge1+inter1+basal2+side2+edge2;

%after
vector_side=[r_a(1,jj-1)-r_b(1,jj), r_a(2,jj-1)-r_b(2,jj), 0];
vector_basal=[r_b(1,jj)-r_b(1,jj-2), r_b(2,jj)-r_b(2,jj-2), 0];
vector_edge=[r_b(1,jj)-virtual(1), r_b(2,jj)-virtual(2), 0];
vector_inter=[r_a(1,jj)-virtual(1), r_a(2,jj)-virtual(2), 0];

basal1=k(2)/2*(norm(vector_basal)-ll_b(jj-2))^2;%+k(2)/2*(sqrt((r_b(1,N)-r_b(1,N-1))^2+(r_b(2,N)-r_b(2,N-1))^2)-ll_b(N-1))^2;
side1=k(3)/2*(norm(vector_side)-ll_s(jj-1))^2;
edge1=k(2)/2*(norm(vector_edge)-ll_b(jj-1))^2;
inter1=k(2)/2*(norm(vector_inter)-ll_s(jj))^2;

vector_side=[r_a(1,jj+1)-r_b(1,jj), r_a(2,jj+1)-r_b(2,jj), 0];
vector_basal=[r_b(1,jj)-r_b(1,jj+2), r_b(2,jj)-r_b(2,jj+2), 0];
vector_edge=[r_b(1,jj)-virtual(1), r_b(2,jj)-virtual(2), 0];

basal2=k(2)/2*(norm(vector_basal)-ll_b(jj+1))^2;%+k(2)/2*(sqrt((r_b(1,N)-r_b(1,N-1))^2+(r_b(2,N)-r_b(2,N-1))^2)-ll_b(N-1))^2;
side2=k(3)/2*(norm(vector_side)-ll_s(jj+1))^2;
edge2=k(2)/2*(norm(vector_edge)-ll_b(jj))^2;

after=basal1+side1+edge1+inter1+basal2+side2+edge2;

spring=after-before;

end

if ID_edge(jj)==0.75 || ID_edge(jj)==1.75
    
    virtual(1)=abs((r_b(1,jj+1)-r_b(1,jj))/2+r_b(1,jj));
 virtual(2)=abs((r_b(2,jj+1)-r_b(2,jj))/2+r_b(2,jj));
 virtual_real=(theta(jj+1)-theta(jj))/2+theta(jj);
    
%before
vector_side=[r_a(1,jj-1)-r_b(1,jj-1), r_a(2,jj-1)-r_b(2,jj-1), 0];
vector_basal=[r_b(1,jj-1)-r_b(1,jj-2), r_b(2,jj-1)-r_b(2,jj-2), 0];
vector_edge=[r_b(1,jj)-r_b(1,jj-1), r_b(2,jj)-r_b(2,jj-1), 0];
vector_inter=[r_b(1,jj+1)-r_b(1,jj), r_b(2,jj+1)-r_b(2,jj), 0];

basal1=k(2)/2*(norm(vector_basal)-ll_b(jj-2))^2;%+k(2)/2*(sqrt((r_b(1,N)-r_b(1,N-1))^2+(r_b(2,N)-r_b(2,N-1))^2)-ll_b(N-1))^2;
side1=k(3)/2*(norm(vector_side)-ll_s(jj-1))^2;
edge1=k(2)/2*(norm(vector_edge)-ll_b(jj-1))^2;
inter1=k(2)/2*(norm(vector_inter)-ll_b(jj))^2;


before=basal1+side1+edge1+inter1;

%after
vector_side=[r_a(1,jj-1)-r_b(1,jj), r_a(2,jj-1)-r_b(2,jj), 0];
vector_basal=[r_b(1,jj)-r_b(1,jj-2), r_b(2,jj)-r_b(2,jj-2), 0];
vector_edge=[r_b(1,jj)-virtual(1), r_b(2,jj)-virtual(2), 0];
vector_inter=[r_b(1,jj+1)-virtual(1), r_b(2,jj+1)-virtual(2), 0];

basal1=k(2)/2*(norm(vector_basal)-ll_b(jj-2))^2;%+k(2)/2*(sqrt((r_b(1,N)-r_b(1,N-1))^2+(r_b(2,N)-r_b(2,N-1))^2)-ll_b(N-1))^2;
side1=k(3)/2*(norm(vector_side)-ll_s(jj-1))^2;
edge1=k(2)/2*(norm(vector_edge)-ll_b(jj-1))^2;
inter1=k(2)/2*(norm(vector_inter)-ll_s(jj))^2;


after=basal1+side1+edge1+inter1;

spring=after-before;
end

if ID_edge(jj)==1.25 || ID_edge(jj)==2.25
    
    virtual(1)=abs((r_b(1,jj-1)-r_b(1,jj))/2+r_b(1,jj));
 virtual(2)=abs((r_b(2,jj-1)-r_b(2,jj))/2+r_b(2,jj));
 
 virtual_real=(theta(jj-1)-theta(jj))/2+theta(jj);
    
%before

vector_inter=[r_b(1,jj-1)-r_b(1,jj), r_b(2,jj-1)-r_b(2,jj), 0];


inter1=k(2)/2*(norm(vector_inter)-ll_b(jj-1))^2;

vector_side=[r_a(1,jj+1)-r_b(1,jj+1), r_a(2,jj+1)-r_b(2,jj+1), 0];
vector_basal=[r_b(1,jj+1)-r_b(1,jj+2), r_b(2,jj+1)-r_b(2,jj+2), 0];
vector_edge=[r_b(1,jj)-r_b(1,jj+1), r_b(2,jj)-r_b(2,jj+1), 0];

basal2=k(2)/2*(norm(vector_basal)-ll_b(jj+1))^2;%+k(2)/2*(sqrt((r_b(1,N)-r_b(1,N-1))^2+(r_b(2,N)-r_b(2,N-1))^2)-ll_b(N-1))^2;
side2=k(3)/2*(norm(vector_side)-ll_s(jj+1))^2;
edge2=k(2)/2*(norm(vector_edge)-ll_b(jj))^2;

before=inter1+basal2+side2+edge2;

%after

vector_inter=[r_b(1,jj-1)-virtual(1), r_b(2,jj-1)-virtual(2), 0];


inter1=k(2)/2*(norm(vector_inter)-ll_b(jj-1))^2;

vector_side=[r_a(1,jj+1)-r_b(1,jj), r_a(2,jj+1)-r_b(2,jj), 0];
vector_basal=[r_b(1,jj)-r_b(1,jj+2), r_b(2,jj)-r_b(2,jj+2), 0];
vector_edge=[r_b(1,jj)-virtual(1), r_b(2,jj)-virtual(2), 0];

basal2=k(2)/2*(norm(vector_basal)-ll_b(jj+1))^2;%+k(2)/2*(sqrt((r_b(1,N)-r_b(1,N-1))^2+(r_b(2,N)-r_b(2,N-1))^2)-ll_b(N-1))^2;
side2=k(3)/2*(norm(vector_side)-ll_s(jj+1))^2;
edge2=k(2)/2*(norm(vector_edge)-ll_b(jj))^2;

after=inter1+basal2+side2+edge2;

spring=after-before;   
    
    
end

energy=+interface+spring;%-grad^2*dtrr
if energy<0
do_it=1;
else
do_it=0;  
end
end

function [delta_L_l, delta_L_r, ID_min_l_e, ID_min_l_n, ID_min_r_e,ID_min_r_n]=interface_length_func(L,ID_initial)
maxval=max(L(ID_initial==1)); %left boundary lengths
        ID_max_l=intersect(find(L==maxval),find(ID_initial==1));%left boundary max
        
        maxval=max(L(ID_initial==2)); %right boundary lengths
        ID_max_r=intersect(find(L==maxval),find(ID_initial==2));%right boundary max
        
       
        minval=min(L(intersect(1:ID_max_l-1,find(ID_initial==1))));
        if isempty(minval)
            ID_min_l_e=length(L);
        else
        ID_min_l_e=intersect(find(L==minval),find(ID_initial==1));
        end
        
        
        minval=min(L(intersect(ID_max_l+1:ID_max_r-1,find(ID_initial==1))));
        if isempty(minval)
            ID_min_l_n=length(L);
        else
        ID_min_l_n=intersect(find(L==minval),find(ID_initial==1));
        end
        
        
        minval=min(L(intersect(ID_max_r+1:length(L),find(ID_initial==2))));
        if isempty(minval)
            ID_min_r_e=length(L);
            
        else
        ID_min_r_e=intersect(find(L==minval),find(ID_initial==2));
        
        end
        
        
        minval=min(L(intersect(ID_max_l+1:ID_max_r-1,find(ID_initial==2))));
        if isempty(minval)
            ID_min_r_n=length(L);
        else
        ID_min_r_n=intersect(find(L==minval),find(ID_initial==2));
        end
        
        
        delta_L_l=-L(ID_min_l_e(1))+L(ID_min_l_n(1));
        delta_L_r=-L(ID_min_r_e(1))+L(ID_min_r_n(1));
        
end
function value= dH_extra(x,y,R,x_0,delta_L,dim,jj,ID_e,ID_n,alpha)
value=0;
if jj==ID_e
if dim==1
    value=alpha*(delta_L)*(-y)*R/(y^2+(x-x_0)^2);
end
if dim==2
    value=alpha*(delta_L)*R/(y^2/(x-x_0)+(x-x_0));
end
end

if jj==ID_n
if dim==1
    value=-alpha*(delta_L)*(-y)*R/(y^2+(x-x_0)^2);
end
if dim==2
    value=-alpha*(delta_L)*R/(y^2/(x-x_0)+(x-x_0));
end
end

end
%right
function value = dx_dr(x,z,x_0,R,z_0,L,theta)

value=cos(theta);

end

function value = dz_dr(x,z,x_0,R,z_0,L,theta)

value=sin(theta);

end

function value=dx_dtheta(R,theta)

value=-R*sin(theta);

end

function value=dz_dtheta(R,theta)

value=R*cos(theta);

end


function [H_tot, H, apical, basal, side, area1, miosin]=Hamiltonian(ll_a,ll_b,ll_s,r,kappa,rr,N,k,gamma,lambda,g,A_0,step,L,radius,x_0_r,x_0_l,delta_L_l,delta_L_r,alpha, mu, ID_initial,Area)
r_a=r(:,1:length(r)/2);
r_b=r(:,1+length(r)/2:end);
for jj=1:length(r_a)-1
vector_side=[r_a(1,jj)-r_b(1,jj), r_a(2,jj)-r_b(2,jj), 0];
vector_apical=[r_a(1,jj+1)-r_a(1,jj), r_a(2,jj+1)-r_a(2,jj), 0];
vector_basal=[r_b(1,jj+1)-r_b(1,jj), r_b(2,jj+1)-r_b(2,jj), 0];
vector_side_2=[r_a(1,jj+1)-r_b(1,jj+1), r_a(2,jj+1)-r_b(2,jj+1), 0];

apical(jj)=k(1)/2*(norm(vector_apical)-ll_a(jj))^2;%+k(1)/2*(sqrt((r_a(1,N)-r_a(1,N-1))^2+(r_a(2,N)-r_a(2,N-1))^2)-ll_a(N-1))^2;
basal(jj)=k(2)/2*(norm(vector_basal)-ll_b(jj))^2;%+k(2)/2*(sqrt((r_b(1,N)-r_b(1,N-1))^2+(r_b(2,N)-r_b(2,N-1))^2)-ll_b(N-1))^2;
side(jj)=k(3)/2*(norm(vector_side)-ll_s(jj))^2;
%attachment=kappa/2*((rr(jj)-1/2*(r_b(1,jj)+r_b(1,jj+1)))^2+(rr(N-1)-1/2*(r_b(1,N)+r_b(1,N-1)))^2);
attachment=0;%-gamma*exp(-1/2*abs(r_b(2,jj)+r_b(2,jj+1))/lambda);%-gamma*exp(-1/2*abs(r_b(2,N)+r_b(2,N-1))/lambda);

sin1=norm(cross(vector_side,vector_basal))/(norm(vector_side)*norm(vector_basal));
sin2=norm(cross(vector_side_2,vector_apical))/(norm(vector_side_2)*norm(vector_apical));
%area1=g/2*abs(0.5*sin1*norm(vector_side)*norm(vector_basal)+0.5*sin2*norm(vector_side_2)*norm(vector_apical)-A_0(jj))^2;
area1(jj)=g/2*(Area(jj)-A_0(jj))^2;
% x=[r_a(1,jj) r_a(1,jj+1) r_b(1,jj) r_b(1,jj+1)];
% y=[r_a(2,jj) r_a(2,jj+1) r_b(2,jj) r_b(2,jj+1)];
% area1=g*(polyarea(x,y)-A_0(jj))^2;

miosin(jj)=norm(vector_apical)*step(jj);
%area1=0;
H(jj)=apical(jj)+basal(jj)+side(jj)+attachment+area1(jj)+miosin(jj);


end
H(length(r_a))=0;
%extra=alpha/2*delta_L_l^2+alpha/2*delta_L_r^2;
extra=-(length(find(ID_initial>0))-1)*mu;
H_tot=sum(H)+extra;
end

function [dH_a, side, front, back, area1, area2, miosin]=dH_a_func(ll_a,ll_s,r_a,r_b,kappa,jj,dim,k,g,delta_x,step,A_0,Area)
if dim ==1
    notdim=2;
else
    notdim=1;
end

if jj==1
    
    side=k(3)*(r_a(dim,jj)-r_b(dim,jj))*(sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2)-ll_s(jj))/sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2);
    front=-k(1)*(r_a(dim,jj+1)-r_a(dim,jj))*(sqrt((r_a(1,jj+1)-r_a(1,jj))^2+(r_a(2,jj+1)-r_a(2,jj))^2)-ll_a(jj))/sqrt((r_a(1,jj+1)-r_a(1,jj))^2+(r_a(2,jj+1)-r_a(2,jj))^2);
    back=0;
    miosin=-(r_a(dim,jj+1)-r_a(dim,jj))*step(jj)/sqrt((r_a(1,jj+1)-r_a(1,jj))^2+(r_a(2,jj+1)-r_a(2,jj))^2);
    area1=1/4*g*(2*r_b(notdim,jj)-2*r_a(notdim,jj+1))*(-1)^(notdim)*(Area(jj)-A_0(jj));
    area2=0;
end




if jj~=1 && jj~=length(r_a)
    side=k(3)*(r_a(dim,jj)-r_b(dim,jj))*(sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2)-ll_s(jj))/sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2);
    back=k(1)*(r_a(dim,jj)-r_a(dim,jj-1))*(sqrt((r_a(1,jj)-r_a(1,jj-1))^2+(r_a(2,jj)-r_a(2,jj-1))^2)-ll_a(jj-1))/sqrt((r_a(1,jj)-r_a(1,jj-1))^2+(r_a(2,jj)-r_a(2,jj-1))^2);  
    front=-k(1)*(r_a(dim,jj+1)-r_a(dim,jj))*(sqrt((r_a(1,jj+1)-r_a(1,jj))^2+(r_a(2,jj+1)-r_a(2,jj))^2)-ll_a(jj))/sqrt((r_a(1,jj+1)-r_a(1,jj))^2+(r_a(2,jj+1)-r_a(2,jj))^2);
    miosin=-(r_a(dim,jj+1)-r_a(dim,jj))*step(jj)/sqrt((r_a(1,jj+1)-r_a(1,jj))^2+(r_a(2,jj+1)-r_a(2,jj))^2)+(r_a(dim,jj)-r_a(dim,jj-1))*step(jj-1)/sqrt((r_a(1,jj)-r_a(1,jj-1))^2+(r_a(2,jj)-r_a(2,jj-1))^2);
    area1=1/4*g*(2*r_b(notdim,jj)-2*r_a(notdim,jj+1))*(-1)^(notdim)*(Area(jj)-A_0(jj));
    area2=1/4*g*(-2*r_b(notdim,jj)+2*r_a(notdim,jj-1))*(-1)^(notdim)*(Area(jj-1)-A_0(jj-1));
end
if jj==51
    a=1;
end
if jj==length(r_a)
    
    side=k(3)*(r_a(dim,jj)-r_b(dim,jj))*(sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2)-ll_s(jj))/sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2);
    front=0;
    back=k(1)*(r_a(dim,jj)-r_a(dim,jj-1))*(sqrt((r_a(1,jj)-r_a(1,jj-1))^2+(r_a(2,jj)-r_a(2,jj-1))^2)-ll_a(jj-1))/sqrt((r_a(1,jj)-r_a(1,jj-1))^2+(r_a(2,jj)-r_a(2,jj-1))^2);
    miosin=+(r_a(dim,jj)-r_a(dim,jj-1))*step(jj-1)/sqrt((r_a(1,jj)-r_a(1,jj-1))^2+(r_a(2,jj)-r_a(2,jj-1))^2);
    area1=0;
    area2=1/4*g*(-2*r_b(notdim,jj)+2*r_a(notdim,jj-1))*(-1)^(notdim)*(Area(jj-1)-A_0(jj-1));
end



dH_a=side+front+back+area1+area2+miosin;
end

function [dH_b,side, front, back, area1, area2]=dH_b_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,dim,k,gamma,lambda, g, delta_x,A_0,Area)
if dim ==1
    notdim=2;
else
    notdim=1;
end

if jj==1
    side=-k(3)*(r_a(dim,jj)-r_b(dim,jj))*(sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2)-ll_s(jj))/sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2);
    front=-k(2)*(r_b(dim,jj+1)-r_b(dim,jj))*(sqrt((r_b(1,jj+1)-r_b(1,jj))^2+(r_b(2,jj+1)-r_b(2,jj))^2)-ll_b(jj))/sqrt((r_b(1,jj+1)-r_b(1,jj))^2+(r_b(2,jj+1)-r_b(2,jj))^2);
    back=0;
    area1=1/4*g*(-2*r_a(notdim,jj)+2*r_b(notdim,jj+1))*(-1)^(notdim)*(Area(jj)-A_0(jj));
    area2=0;

end


if jj~=1 && jj~=length(r_b)
    side=-k(3)*(r_a(dim,jj)-r_b(dim,jj))*(sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2)-ll_s(jj))/sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2);
    front=-k(2)*(r_b(dim,jj+1)-r_b(dim,jj))*(sqrt((r_b(1,jj+1)-r_b(1,jj))^2+(r_b(2,jj+1)-r_b(2,jj))^2)-ll_b(jj))/sqrt((r_b(1,jj+1)-r_b(1,jj))^2+(r_b(2,jj+1)-r_b(2,jj))^2);
    back=k(2)*(r_b(dim,jj)-r_b(dim,jj-1))*(sqrt((r_b(1,jj)-r_b(1,jj-1))^2+(r_b(2,jj)-r_b(2,jj-1))^2)-ll_b(jj-1))/sqrt((r_b(1,jj)-r_b(1,jj-1))^2+(r_b(2,jj)-r_b(2,jj-1))^2); 

    area1=1/4*g*(-2*r_a(notdim,jj)+2*r_b(notdim,jj+1))*(-1)^(notdim)*(Area(jj)-A_0(jj));
    area2=1/4*g*(+2*r_a(notdim,jj)-2*r_b(notdim,jj-1))*(-1)^(notdim)*(Area(jj-1)-A_0(jj-1));
end


if jj==length(r_a)
    side=-k(3)*(r_a(dim,jj)-r_b(dim,jj))*(sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2)-ll_s(jj))/sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2);
    front=0;
    back=k(2)*(r_b(dim,jj)-r_b(dim,jj-1))*(sqrt((r_b(1,jj)-r_b(1,jj-1))^2+(r_b(2,jj)-r_b(2,jj-1))^2)-ll_b(jj-1))/sqrt((r_b(1,jj)-r_b(1,jj-1))^2+(r_b(2,jj)-r_b(2,jj-1))^2);
    area1=0;
    area2=1/4*g*(+2*r_a(notdim,jj)-2*r_b(notdim,jj-1))*(-1)^(notdim)*(Area(jj-1)-A_0(jj-1));
end



if dim==2 && jj==1
    bottom1=gamma/lambda*exp(-1/2*abs(r_b(2,jj)+r_b(2,jj+1))/lambda);
end

if dim==2 && jj==length(r_a)
    bottom2=gamma/lambda*exp(-1/2*abs(r_b(2,jj)+r_b(2,jj-1))/lambda);
    
    %bottom1=-kappa/2*(rr(jj)-1/2*(r_b(1,jj)+r_b(1,jj+1)));
    %bottom2=-kappa/2*(rr(jj-1)-1/2*(r_b(1,jj)+r_b(1,jj-1)));
end

if dim==2 && jj~=1 && jj~=length(r_a)
    bottom1=gamma/lambda*exp(-1/2*abs(r_b(2,jj)+r_b(2,jj+1))/lambda);
    bottom2=gamma/lambda*exp(-1/2*abs(r_b(2,jj)+r_b(2,jj-1))/lambda);
    %bottom1=-kappa/2*(rr(jj)-1/2*(r_b(1,jj)+r_b(1,jj+1)));
    %bottom2=-kappa/2*(rr(jj-1)-1/2*(r_b(1,jj)+r_b(1,jj-1)));
end
bottom1=0;
bottom2=0;
dH_b=side+front+back+bottom1+bottom2+area1+area2;


end

function dp = divpro(r,lmin,lmax,dt,gamma0)

    % probability for a cell to divide
    % zero below lmin, constant above lmax, linear between lmin and lmax
    % output:   0 cell does not divide in dt
    %           1 cell divides
    ll = r(2:end)-r(1:end-1);
    N = length(ll);
    gamma = (ll-lmin)/(lmax-lmin);
    gamma(ll<lmin) = 0;
    gamma(ll>lmin) = 1;
    %gamma(l>lmax) = 1;
    
    gamma = gamma*dt*gamma0;
    
    p = rand(N,1);
    dp = zeros(N,1);
    %dp(gamma-p>0) = 1; % setting cell division to 0

end


function [rp,rrp,lp,IDp] = division(r,rr,ll,ID,CellDivides)

    ix = find(CellDivides == 1)';
    ix = ix + (0:1:length(ix)-1); %shift index because will add cells
    for i=ix
        N= length(r);
        rp(1:i,1)     = r(1:i);
        rp(i+1,1)     = 0.5*(r(i)+r(i+1));
        rp(i+2:N+1,1) = r(i+1:N);
 
        rrp(1:i,1)     = rr(1:i);
        rrp(i,  1)     = 0.5*(rp(i)+rp(i+1));
        rrp(i+1,1)     = 0.5*(rp(i+1)+rp(i+2));
        if (i < N-1)
            rrp(i+2:N,1) = rr(i+1:N-1);
        end
        
        IDp(1:i,1) = ID(1:i);
        IDp(i,1)   = ID(i); % new cells get same ecad\ncad classification as mother
        IDp(i+1,1) = ID(i);
        if (i < N-1)
            IDp(i+2:N,1) = ID(i+1:N-1);
        end
        
        lp(1:i,1) = ll(1:i);
        lp(i,1)   = ll(i)/2; % new cells get half the area
        lp(i+1,1) = ll(i)/2;
        if (i < N-1)
            lp(i+2:N,1) = ll(i+1:N-1);
        end
        r = rp;
        rr= rrp;
        ll = lp;
        ID=IDp;
    end
end

function dhl_a=dhlfunc_a(ll_a,kappa,r_b,fext) % insert cell length, attachement strength, attachment point and force


    dhl_a = zeros(length(ll_a)+1,1);    
    dhl_a(1:end) = r_b;
    %dhl(2:end)   = dhl(2:end) + kappa*0.5*rr;
    dhl_a(1:end-1) = dhl_a(1:end-1) - ll_a;
    dhl_a(2:end)   = dhl_a(2:end)   + ll_a;
    
    %external forces
    dhl_a(1)   = dhl_a(1)   - fext;
    dhl_a(end) = dhl_a(end) + fext;
    
end

function dhl_b=dhlfunc_b(ll_b,kappa,rr,fext,r_a) % insert cell length, attachement strength, attachment point and force


    dhl_b = zeros(length(ll_b)+1,1);
    dhl_b(1:end)=r_a;
    dhl_b(1:end-1) =dhl_b(1:end-1)+ kappa*0.5*rr;
    dhl_b(2:end)   = dhl_b(2:end) + kappa*0.5*rr;
    dhl_b(1:end-1) = dhl_b(1:end-1) - ll_b;
    dhl_b(2:end)   = dhl_b(2:end)   + ll_b;
    
    %external forces
    dhl_b(1)   = dhl_b(1)   - fext;
    dhl_b(end) = dhl_b(end) + fext;
    
end

function [dhs_a] = dhfunc_a(r,kappa) % insert edges and attachement strength


    N = length(r);

    % inner part
    i = 2:N-1;
    j = 2:N-1;
    s = ones(1,N-2)*( 3); 
    
    i = [i (2:N-1)-1];
    j = [j (2:N-1)];
    s = [s ones(1,N-2)*(-1 )];

    i = [i (2:N-1)+1];
    j = [j (2:N-1)];
    s = [s ones(1,N-2)*(-1 )];
    
    % boundary
    i = [i 1];
    j = [j 1];
    s = [s  2];

    i = [i 2];
    j = [j 1];
    s = [s -1 ];

    i = [i N];
    j = [j N];
    s = [s  2];

    i = [i N-1];
    j = [j N];
    s = [s -1 ];    
    dhs_a = sparse(i,j,s,N,N);    
    
end

function [dhs_b] = dhfunc_b(r,kappa) % insert edges and attachement strength


    N = length(r);

    % inner part
    i = 2:N-1;
    j = 2:N-1;
    s = ones(1,N-2)*( 3 + kappa/2); 
    
    i = [i (2:N-1)-1];
    j = [j (2:N-1)];
    s = [s ones(1,N-2)*(-1 + kappa/4)];

    i = [i (2:N-1)+1];
    j = [j (2:N-1)];
    s = [s ones(1,N-2)*(-1 + kappa/4)];
    
    % boundary
    i = [i 1];
    j = [j 1];
    s = [s  2 + kappa/4];

    i = [i 2];
    j = [j 1];
    s = [s -1 + kappa/4];

    i = [i N];
    j = [j N];
    s = [s  2 + kappa/4];

    i = [i N-1];
    j = [j N];
    s = [s -1 + kappa/4];    
    dhs_b = sparse(i,j,s,N,N);    
    
end

function [r_a,r_b,rr,ll_a,ll_b,ll_s,ID,A_0,step,ID_initial,ID_edge, ID_v] = ChainModelInitial(N, sigma,lmax,m,ll_b_e_g,ll_b_n_g,Nnumber) % Number of cells and noise in initial condition
    
    rng(0);
    r_a(1,:) = ((1+N-1:1:N+N-1)'-(N+1)/2);%  + rand(N,1)*sigma;% define cell center, initially at the 0.5, for example one cell would have center 0.5, then add noise
    
    rng(1)
    r_a(2,:)=ones(N,1);%+rand(N,1)*sigma;
 
    rng(2);
    r_b(1,:) = ((1+N-1:1:N+N-1)'-(N+1)/2);%  + rand(N,1)*sigma; % define cell center, initially at the 0.5, for example one cell would have center 0.5, then add noise
    
    
    rng(3)
    r_b(2,:)=zeros(N,1);%rand(N,1)*sigma
    
    % ID_initial marks the apical point of the boundary between N/E
    ID_initial=zeros(length(r_a),1);
    ID_initial((N-1-Nnumber)/2+1)= 1;
    ID_initial((N-1-Nnumber)/2+Nnumber+1)= 2;
    ID_initial(1)= -1;
    ID_initial(end)= -1;
    rr = 0.5*(r_b(1,2:end)+r_b(1,1:end-1));% + rand(N-1,1)*sigma; %attachment points (on the substrate)
    %l  = r(2:end)-r(1:end-1);%ones(N-1,1);
    %l = l-min(l)+.2;
    ll_a=1*ones(N-1,1); %ll_a = unifrnd(1,lmax,N-1,1)*4; %random numbers from distribution, to get size of cells %ecad top
    %ll_a(length(ll_a)/4+1:length(ll_a)*3/4)=1*ones(Nnumber,1);%ll_a(length(ll_a)/4+1:length(ll_a)*3/4)=unifrnd(1,lmax,(N-1)/2,1); %ncad top 0.1 0.5
    
    ll_b=ll_b_e_g*ones(N-1,1);%ll_b = unifrnd(1,lmax,N-1,1); %random numbers from distribution, to get size of cells %ecad bottom
    %ll_b(length(ll_b)/4+1:length(ll_b)*3/4)=ll_b_n_g*ones(Nnumber,1);%ll_b(length(ll_b)/4+1:length(ll_b)*3/4)=unifrnd(1,lmax,(N-1)/2,1); %ncad bottom
    
    ll_s=1*ones(N,1);%ll_b = unifrnd(1,lmax,N-1,1); %random numbers from distribution, to get size of cells %ecad bottom
    %ll_s((length(ll_s)-1)/4+1:(length(ll_s)-1)*3/4+1)=1*ones((N-1)/2+1,1);
    ID=ones(length(ll_a),1);
    ID((N-1-Nnumber)/2+1:(N-1-Nnumber)/2+Nnumber)= 2;
    
    A_0=ones(length(ll_a));
    
    step=0*ones(N-1,1);
    step((N-1-Nnumber)/2+1:(N-1-Nnumber)/2+Nnumber)=m*ones(Nnumber,1);
    
    
    ID_edge=zeros(length(r_a),1);
    ID_edge((N-1-Nnumber)/2+1)= 1;
    ID_edge((N-1-Nnumber)/2+Nnumber+1)= 2;
    
       ID_v=ones(length(r_a),1);
    ID_v((N-1-Nnumber)/2+1:(N-1-Nnumber)/2+Nnumber+1)= 2;
   
end