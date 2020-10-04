function [points,cells] = PlasticGrowthH2D_line(K,gamma,g,m,alpha)

    
    
    % spring mechancis model for tissue/colony growth
    % r     positions of cells
    % rr    attachment points
    % ll     length of a cell
    % kappa attachment strength

    
    PlasticGrowth = 0; % Growth in response to stress = 1; Constant Growth = 0;
    
    % parameters:
            
       
        NTimes = 400000; % Number of simulation runs;
        Plots  = 1000; % When to plot distribution
        
        Ncells = 41;   % number of edges at start
        sigma  = 1;    % standard deviation of initial condition
        kappa  = 1;    % attachment strength
        gamma0 = 0.04;   % cell division rate 
        eta    = 0;    % cell growth rate %0.08 %0 is set for no cell growth
        nu     = 1;   % friction coefficient 
        dtrr   = 0.001; % time step for attachment point relaxation
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
        radius=10;
        x_0=30.5;
        x_offset=0;
        alpha=alpha;

    % Chain initiation, now we also define the line in the N/E boundary
    [r_a,r_b,rr,ll_a,ll_b,ll_s,ID,A_0,step,ID_initial,ID_edge] = ChainModelInitial(Ncells,sigma,lmax,m); % from this we get cell edges, attachement points and cell length
    r_a(1,11)=x_0-sqrt(radius^2-1);
    r_a(1,31)=sqrt(radius^2-1)+x_0;
    k=[k_a k_b k_s];
    
    H=Hamiltonian(ll_a,ll_b,ll_s,[r_a,r_b],kappa,rr,Ncells,k,gamma,lambda,g,A_0,step,0,0,0,0,0,0,0);
    
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
    points(1).radius=10;
    points(1).x_0_r=0;
    points(1).x_0_l=0;
    cells = zeros(1,NTimes);
    %divposdt = [];
    cmp = jet(NTimes/Plots+1);
    idx = 1;
    
   
%     dhs_a = dhfunc_a(r_a, kappa);
%     dhs_b = dhfunc_b(r_b, kappa);
    L=zeros(length(r_a),1);
    L_a=zeros(length(r_a),1);
    L_a(ID_edge==1)=  -radius*atan(r_a(2,find(ID_edge==1))/(r_a(1,find(ID_edge==1))-x_0));
    L_a(ID_edge==2)= radius*atan(r_a(2,find(ID_edge==2))/(r_a(1,find(ID_edge==2))-x_0));
    for i= 1:NTimes
       
        n=0;
        
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
        
        for jj=1:Ncells
            
            n=n+1;
            if jj~=length(r_a)
            x1=[r_b(1,jj) r_b(1,jj+1) r_a(1,jj+1) r_a(1,jj)];
            y1=[r_b(2,jj) r_b(2,jj+1) r_a(2,jj+1) r_a(2,jj)];
            Area(jj)=polyarea(x1,y1);
            end
            
            dH_a_x(1,n)=dH_a_func(ll_a,ll_s,r_a,r_b,kappa,jj,1,k,g,delta_x,step,A_0,Area);
            dH_a_z(1,n)=dH_a_func(ll_a,ll_s,r_a,r_b,kappa,jj,2,k,g,delta_x,step,A_0,Area);
            
            dH_b_x(1,n)=dH_b_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,1,k,gamma,lambda,g,delta_x,A_0,Area);
            dH_b_z(1,n)=0;%dH_b_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,2,k,gamma,lambda,g,delta_x,A_0,Area);  
            
%           [dA, A2]=dArea_func(dH_a_x(1,n),dH_a_z(1,n),dH_b_x(1,n),dH_b_z(1,n),r_a,r_b,jj,g,dtrr);
%           Area(i,jj)=A2;
%           dArea(n,1:4)=dA;
            dH_l(1,n)=0;
            dH_r(1,n)=0;
            dH_x_0(1,n)=0;
             if ID_initial(jj)==1
                dH_b_z(1,n)=dH_b_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,2,k,gamma,lambda,g,delta_x,A_0,Area)+dH_extra(r_b(1,jj),r_b(2,jj),radius,x_0,delta_L_l,2,jj,ID_min_l_e,ID_min_l_n,alpha);
                dH_b_x(1,n)=dH_b_x(1,n)+dH_extra(r_b(1,jj),r_b(2,jj),radius,x_0,delta_L_l,1,jj,ID_min_l_e,ID_min_l_n,alpha);
                
                dx=-dx_dl(r_b(1,jj),r_b(2,jj),x_0,radius);
                dz=dz_dl(r_b(1,jj),r_b(2,jj),x_0,radius);
                dH_l(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
            
                dx=-dx_dr(r_b(1,jj),r_b(2,jj),x_0,radius);
                dz=dz_dr(r_b(1,jj),r_b(2,jj),x_0,radius);
                dH_r(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
                
                dx=-1;
                dz=0;
                dH_x_0(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
                if ID_edge(jj)==1
%                 dH_a_z(1,n)=dH_a_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,2,k,gamma,lambda,g,delta_x,A_0,Area);
%                 dH_a_x(1,n)=dH_a_x(1,n);   
                    dx=-dx_dl(r_a(1,jj),r_a(2,jj),x_0,radius);
                    dz=dz_dl(r_a(1,jj),r_a(2,jj),x_0,radius);
                    dH_l_a(1,n)=dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                    
                    dx=-dx_dr(r_a(1,jj),r_a(2,jj),x_0,radius);
                    dz=dz_dr(r_a(1,jj),r_a(2,jj),x_0,radius);
                    dH_r(1,n)=dH_r(1,n)+ dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                
                    dx=-1;
                    dz=0;
                    dH_x_0(1,n)= dH_x_0(1,n)+ dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                
                
                end
                
            end
            if ID_initial(jj)==2
                dH_b_z(1,n)=dH_b_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,2,k,gamma,lambda,g,delta_x,A_0,Area)+dH_extra(r_b(1,jj),r_b(2,jj),radius,x_0,delta_L_r,2,jj,ID_min_r_e,ID_min_r_n,alpha); 
                dH_b_x(1,n)=dH_b_x(1,n)+dH_extra(r_b(1,jj),r_b(2,jj),radius,x_0,delta_L_r,1,jj,ID_min_r_e,ID_min_r_n,alpha);
                
                dx=dx_dl(r_b(1,jj),r_b(2,jj),x_0,radius);
                dz=dz_dl(r_b(1,jj),r_b(2,jj),x_0,radius);
                dH_l(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
            
                dx=dx_dr(r_b(1,jj),r_b(2,jj),x_0,radius);
                dz=dz_dr(r_b(1,jj),r_b(2,jj),x_0,radius);
                dH_r(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
                
                dx=1;
                dz=0;
                dH_x_0(1,n)= dH_b_x(1,jj)*dx + dH_b_z(1,jj)*dz;
                if ID_edge(jj)==2
%                 dH_a_z(1,n)=dH_a_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,2,k,gamma,lambda,g,delta_x,A_0,Area);
%                 dH_a_x(1,n)=dH_a_x(1,n);   
                    dx=dx_dl(r_a(1,jj),r_a(2,jj),x_0,radius);
                    dz=dz_dl(r_a(1,jj),r_a(2,jj),x_0,radius);
                    dH_l_a(1,n)=dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                    
                    dx=dx_dr(r_a(1,jj),r_a(2,jj),x_0,radius);
                    dz=dz_dr(r_a(1,jj),r_a(2,jj),x_0,radius);
                    dH_r(1,n)=dH_r(1,n)+ dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                
                    dx=1;
                    dz=0;
                    dH_x_0(1,n)= dH_x_0(1,n)+ dH_a_x(1,jj)*dx + dH_a_z(1,jj)*dz;
                
                
                end
            end
        end
        dH_r_total=sum(dH_r);
        dH_x_0_total=sum(dH_x_0);


        N=Ncells;
        cells(i) = Ncells;
        fullgrad=[dH_a_x(1,find(ID_edge==0)) dH_a_z(1,find(ID_edge==0)) dH_b_x(1,find(ID_initial==0)) dH_b_z(1,find(ID_initial==0)) dH_l dH_r_total dH_x_0_total dH_l_a];
        gradnorm=norm(fullgrad);
        %rest
        
        radius=radius+dtrr/10*(-dH_r_total)/gradnorm;
        x_offset=x_offset+dtrr/10*(-dH_x_0_total)/gradnorm;
        x_0_r=x_0+x_offset;
        x_0_l=x_0-x_offset;
        
        for jj=2:length(r_a)-1
           
            
            
            if ID_initial(jj)==2

                L(jj)=L(jj)+dtrr*(-dH_l(1,jj))/gradnorm;
                if L(jj)<0
                    L(jj)=0;
                end
                r_b(1,jj)=radius*cos(L(jj)/radius) + x_0_r;
                r_b(2,jj)=radius*sin(L(jj)/radius);
                r_b(2,(r_b(2,:)<0))=0;
                if ID_edge(jj)==2
                L_a(jj)=L_a(jj)+dtrr*(-dH_l_a(1,jj))/gradnorm;
                
                r_a(1,jj)=radius*cos(L_a(jj)/radius) + x_0_r;
                r_a(2,jj)=radius*sin(L_a(jj)/radius);   
   
                else
                r_a(1,jj)=r_a(1,jj)+dtrr*(-dH_a_x(1,jj))/gradnorm;%+dArea(:,1)');%+dtrr*noise
                r_a(2,jj)=r_a(2,jj)+dtrr*(-dH_a_z(1,jj))/gradnorm;
                end
            end
            
            if ID_initial(jj)==1
               
                
                L(jj)=L(jj)+dtrr*(-dH_l(1,jj))/gradnorm;
                if L(jj)<0
                    L(jj)=0;
                end
                
                r_b(1,jj)= x_0_l-radius*cos(L(jj)/radius);
                r_b(2,jj)=radius*sin(L(jj)/radius);
                r_b(2,(r_b(2,:)<0))=0;
                if ID_edge(jj)==1
                L_a(jj)=L_a(jj)+dtrr*(-dH_l_a(1,jj))/gradnorm;
                
                r_a(1,jj)=x_0_l-radius*cos(L_a(jj)/radius);
                r_a(2,jj)=radius*sin(L_a(jj)/radius);   
   
                else
                r_a(1,jj)=r_a(1,jj)+dtrr*(-dH_a_x(1,jj))/gradnorm;%+dArea(:,1)');%+dtrr*noise
                r_a(2,jj)=r_a(2,jj)+dtrr*(-dH_a_z(1,jj))/gradnorm;
                end
            end
            
            if ID_initial(jj)==0
            noise=.2*randn(length(r_a),1);
            r_a(1,jj)=r_a(1,jj)+dtrr*(-dH_a_x(1,jj))/gradnorm;%+dArea(:,1)');%+dtrr*noise
            r_a(2,jj)=r_a(2,jj)+dtrr*(-dH_a_z(1,jj))/gradnorm;%+dArea(:,2)');
        
            noise=.2*randn(length(r_b),1);
            r_b(1,jj)=r_b(1,jj)+dtrr*(-dH_b_x(1,jj))/gradnorm;%+dArea(:,3)');%+dtrr*noise?
            r_b(2,jj)=r_b(2,jj);%+dtrr*(-dH_b_z(1,jj))/gradnorm;%+dArea(:,4)');
            r_b(2,(r_b(2,:)<0))=0;
            if 0.25>=abs(x_0_l-radius-r_b(1,jj))
                ID_initial(jj)=1;
                 r_b(1,jj)=x_0_l-radius;
            end
            if 0.25>=abs(x_0_r+radius-r_b(1,jj))
                ID_initial(jj)=2;
                r_b(1,jj)=x_0_r+radius;
            end
            end 
        end
        %adhesion relaxation
        xi      = .2*randn(length(rr),1);        %Generate random noise
        xi(1)   =  -abs(xi(1))-fext;        %active traction on boundaries
        xi(end) = abs(xi(end))+fext;
        rr      = rr + dtrr  *(kappa/nu * (0.5*(r_b(1,1:end-1)+r_b(1,2:end))-rr) ); %+ sqrt(dtrr)*xi;
        %rr      = rr + dtrr  * (kappa/nu * (r-rr) ) + sqrt(dtrr)*xi; %spring from basal point to apical point above it
        
        H=Hamiltonian(ll_a,ll_b,ll_s,[r_a,r_b],kappa,rr,Ncells,k,gamma,lambda,g,A_0,step,L,radius,x_0_r,x_0_l,delta_L_l,delta_L_r,alpha);
        
        % Growth of cells:
        if PlasticGrowth == 1   % StressDependant
            
            s     = r_a(2:end)-r_a(1:end-1)-ll_a;
            sp    = 0.;
            ix    = find(s>sp); 
            ll_a(ix) = ll_a(ix) + dtrr*eta; %growth of cell
        else                    % Constant
            
            ll_a     = ll_a + dtrr*eta;
        end
        
        %ll_a(ll_a>lmax) = lmax; % Cells don't grow above lmax.
        if i~=1

            dArea=abs(A_0-Area)./A_0*100;
        end
            
        if mod(i,Plots) == 0
           
            points(i/Plots+1).q_a = r_a ;
            points(i/Plots+1).q_b = r_b ;
            points(i/Plots+1).R = rr;
            points(i/Plots+1).l_a = ll_a;
            points(i/Plots+1).l_b = ll_b;
            points(i/Plots+1).l_s = ll_s;
            points(i/Plots+1).ID=ID;
            points(i/Plots+1).H=H;
            points(i/Plots+1).dH_a=dH_a_x;
            points(i/Plots+1).dH_b=dH_b_x;
            points(i/Plots+1).dH_a_z=dH_a_z;
            points(i/Plots+1).dH_b_z=dH_b_z;
            points(i/Plots+1).dArea=dArea;
            points(i/Plots+1).rmsA=[];
            points(i/Plots+1).radius=radius;
            points(i/Plots+1).x_0_r=x_0_r;
            points(i/Plots+1).x_0_l=x_0_l;
           
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
   
end
%         if i~=1
%             
%            deltaA=(Area(i,:)-Area(i-1,:))./Area(i,:);
%            rmsA=sqrt(nanmean(deltaA(:).^2));
%         
%         else
%         rmsA=[]; 
%         end
%         dhl_a=dhlfunc_a(ll_a,kappa,r_b,fext);
%         r_a_min = dhs_a\dhl_a; % this computes new edges(?)
%         
%         dhl_b=dhlfunc_b(ll_b,kappa,rr,fext,r_a);
%         r_b_min=dhs_b\dhl_b;

% minimization attempt
% fun=@(x)Hamiltonian(ll_a,ll_b,ll_s,x,kappa,rr,Ncells,k,gamma,lambda,g,A_0,step);
% options=optimset('MaxFunEvals', 2000*numel([r_a,r_b]));
% x=fminsearch(fun,[r_a,r_b],options);

% dH_a_x=x(1,1:length(x)/2)-r_a(1,:);
% dH_a_z=x(2,1:length(x)/2)-r_a(2,:);
% dH_b_x=x(1,1+length(x)/2:end)-r_b(1,:);
% dH_b_z=x(2,1+length(x)/2:end)-r_b(2,:);
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
function value = dx_dr(x,z,x_0,R)

value=cos(atan(z/(x-x_0)));

end

function value = dz_dr(x,z,x_0,R)

value=sin(atan(z/(x-x_0)));

end

function value=dx_dl(x,z,x_0,R)

value=-sin(atan(z/(x-x_0)));

end

function value=dz_dl(x,z,x_0,R)

value=cos(atan(z/(x-x_0)));

end
% function [dArea, A2]=dArea_func(dH_a_x,dH_a_z,dH_b_x,dH_b_z,r_a,r_b,jj,g,dtrr)
% 
% if jj==1
% 
% x2=[r_a(1,jj) r_a(1,jj+1) r_b(1,jj+1) r_b(1,jj)];
% y2=[r_a(2,jj) r_a(2,jj+1) r_b(2,jj+1) r_b(2,jj)];
% 
% 
% A2=polyarea(x2,y2);
% 
% delta_x=-dtrr*dH_a_x;
% dx=[delta_x 0 0 0];
% area1=g*(polyarea(x2+dx,y2)-A2);
% 
% delta_x=-dtrr*dH_a_z;
% dx=[delta_x 0 0 0];
% area2=-g*(polyarea(x2,y2+dx)-A2);
% 
% delta_x=-dtrr*dH_b_x;
% dx=[0 0 0 delta_x];
% area3=g*(polyarea(x2+dx,y2)-A2);
% 
% delta_x=-dtrr*dH_b_z;
% dx=[0 0 0 delta_x];
% area4=g*(polyarea(x2,y2+dx)-A2);
% 
% 
% end
% if jj==length(r_a)
% x1=[r_a(1,jj) r_a(1,jj-1) r_b(1,jj-1) r_b(1,jj)];
% y1=[r_a(2,jj) r_a(2,jj-1) r_b(2,jj-1) r_b(2,jj)];
% 
% 
% A1=polyarea(x1,y1);
% A2=NaN;
% 
% delta_x=-dtrr*dH_a_x;
% dx=[delta_x 0 0 0];
% area1=-g*(polyarea(x1+dx,y1)-A1) ;
% 
% delta_x=-dtrr*dH_a_z;
% dx=[delta_x 0 0 0];
% area2=-g*(polyarea(x1,y1+dx)-A1) ;
% 
% delta_x=-dtrr*dH_b_x;
% dx=[0 0 0 delta_x];
% area3=-g*(polyarea(x1+dx,y1)-A1) ;
% 
% delta_x=-dtrr*dH_b_z;
% dx=[0 0 0 delta_x];
% area4=g*(polyarea(x1,y1+dx)-A1) ;        
%         
%     
%     
% end
% if jj~=1 && jj~=length(r_a)
% x1=[r_a(1,jj) r_a(1,jj-1) r_b(1,jj-1) r_b(1,jj)];
% y1=[r_a(2,jj) r_a(2,jj-1) r_b(2,jj-1) r_b(2,jj)];
% 
% x2=[r_a(1,jj) r_a(1,jj+1) r_b(1,jj+1) r_b(1,jj)];
% y2=[r_a(2,jj) r_a(2,jj+1) r_b(2,jj+1) r_b(2,jj)];
% 
% A1=polyarea(x1,y1);
% A2=polyarea(x2,y2);   
%    
% delta_x=-dtrr*dH_a_x;
% dx=[delta_x 0 0 0];
% area1=-g*(polyarea(x1+dx,y1)-A1) +g*(polyarea(x2+dx,y2)-A2);
% 
% delta_x=-dtrr*dH_a_z;
% dx=[delta_x 0 0 0];
% area2=-g*(polyarea(x1,y1+dx)-A1) -g*(polyarea(x2,y2+dx)-A2);
% 
% delta_x=-dtrr*dH_b_x;
% dx=[0 0 0 delta_x];
% area3=-g*(polyarea(x1+dx,y1)-A1) +g*(polyarea(x2+dx,y2)-A2);
% 
% delta_x=-dtrr*dH_b_z;
% dx=[0 0 0 delta_x];
% area4=g*(polyarea(x1,y1+dx)-A1) +g*(polyarea(x2,y2+dx)-A2);    
%     
% end
% 
% 
% 
% dArea=[area1 area2 area3 area4];
% dArea(isnan(dArea))=0;
% end



function H=Hamiltonian(ll_a,ll_b,ll_s,r,kappa,rr,N,k,gamma,lambda,g,A_0,step,L,radius,x_0_r,x_0_l,delta_L_l,delta_L_r,alpha)
r_a=r(:,1:length(r)/2);
r_b=r(:,1+length(r)/2:end);
for jj=1:length(r_a)-1
vector_side=[r_a(1,jj)-r_b(1,jj), r_a(2,jj)-r_b(2,jj), 0];
vector_apical=[r_a(1,jj+1)-r_a(1,jj), r_a(2,jj+1)-r_a(2,jj), 0];
vector_basal=[r_b(1,jj+1)-r_b(1,jj), r_b(2,jj+1)-r_b(2,jj), 0];
vector_side_2=[r_a(1,jj+1)-r_b(1,jj+1), r_a(2,jj+1)-r_b(2,jj+1), 0];

apical=k(1)/2*(norm(vector_apical)-ll_a(jj))^2;%+k(1)/2*(sqrt((r_a(1,N)-r_a(1,N-1))^2+(r_a(2,N)-r_a(2,N-1))^2)-ll_a(N-1))^2;
basal=k(2)/2*(norm(vector_basal)-ll_b(jj))^2;%+k(2)/2*(sqrt((r_b(1,N)-r_b(1,N-1))^2+(r_b(2,N)-r_b(2,N-1))^2)-ll_b(N-1))^2;
side=k(3)/2*(norm(vector_side)-ll_s(jj))^2;
%attachment=kappa/2*((rr(jj)-1/2*(r_b(1,jj)+r_b(1,jj+1)))^2+(rr(N-1)-1/2*(r_b(1,N)+r_b(1,N-1)))^2);
attachment=-gamma*exp(-1/2*abs(r_b(2,jj)+r_b(2,jj+1))/lambda);%-gamma*exp(-1/2*abs(r_b(2,N)+r_b(2,N-1))/lambda);

sin1=norm(cross(vector_side,vector_basal))/(norm(vector_side)*norm(vector_basal));
sin2=norm(cross(vector_side_2,vector_apical))/(norm(vector_side_2)*norm(vector_apical));
area1=g/2*abs(0.5*sin1*norm(vector_side)*norm(vector_basal)+0.5*sin2*norm(vector_side_2)*norm(vector_apical)-A_0(jj))^2;
% x=[r_a(1,jj) r_a(1,jj+1) r_b(1,jj) r_b(1,jj+1)];
% y=[r_a(2,jj) r_a(2,jj+1) r_b(2,jj) r_b(2,jj+1)];
% area1=g*(polyarea(x,y)-A_0(jj))^2;

miosin=norm(vector_apical)*step(jj);
%area1=0;
H(jj)=apical+basal+side+attachment+area1+miosin;


end
extra=alpha/2*delta_L_l^2+alpha/2*delta_L_r^2;
H=sum(H)+extra;
end

% function dH_a=dH_a_func_1(ll_a,ll_s,r_a,r_b,kappa,jj,dim,k)
% if jj==1
% side=k(3)*(r_a(dim,jj)-r_b(dim,jj))*(sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2)-ll_s(jj))/sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2);
% front=-k(1)*(r_a(dim,jj+1)-r_a(dim,jj))*(sqrt((r_a(1,jj+1)-r_a(1,jj))^2+(r_a(2,jj+1)-r_a(2,jj))^2)-ll_a(jj))/sqrt((r_a(1,jj+1)-r_a(1,jj))^2+(r_a(2,jj+1)-r_a(2,jj))^2);
% %back=k(1)*(r_a(dim,jj)-r_a(dim,jj-1))*(sqrt((r_a(1,jj)-r_a(1,jj-1))^2+(r_a(2,jj)-r_a(2,jj-1))^2)-ll_a(jj-1))/sqrt((r_a(1,jj)-r_a(1,jj-1))^2+(r_a(2,jj)-r_a(2,jj-1))^2);
% dH_a=side+front;
% end
% 
% if jj==length(r_a)
% side=k(3)*(r_a(dim,jj)-r_b(dim,jj))*(sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2)-ll_s(jj))/sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2);
% %front=-k(1)*(r_a(dim,jj+1)-r_a(dim,jj))*(sqrt((r_a(1,jj+1)-r_a(1,jj))^2+(r_a(2,jj+1)-r_a(2,jj))^2)-ll_a(jj))/sqrt((r_a(1,jj+1)-r_a(1,jj))^2+(r_a(2,jj+1)-r_a(2,jj))^2);
% back=k(1)*(r_a(dim,jj)-r_a(dim,jj-1))*(sqrt((r_a(1,jj)-r_a(1,jj-1))^2+(r_a(2,jj)-r_a(2,jj-1))^2)-ll_a(jj-1))/sqrt((r_a(1,jj)-r_a(1,jj-1))^2+(r_a(2,jj)-r_a(2,jj-1))^2);
% dH_a=side+back;    
%     
%     
%     
% end
% end

% function dH_b=dH_b_func_1(ll_b,ll_s,r_a,r_b,kappa,jj,rr,dim,k,gamma,lambda)
% if jj==1
% side=-k(3)*(r_a(dim,jj)-r_b(dim,jj))*(sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2)-ll_s(jj))/sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2);
% front=-k(2)*(r_b(dim,jj+1)-r_b(dim,jj))*(sqrt((r_b(1,jj+1)-r_b(1,jj))^2+(r_b(2,jj+1)-r_b(2,jj))^2)-ll_b(jj))/sqrt((r_b(1,jj+1)-r_b(1,jj))^2+(r_b(2,jj+1)-r_b(2,jj))^2);
% %back=k(2)*(r_b(dim,jj)-r_b(dim,jj-1))*(sqrt((r_b(1,jj)-r_b(1,jj-1))^2+(r_b(2,jj)-r_b(2,jj-1))^2)-ll_b(jj-1))/sqrt((r_b(1,jj)-r_b(1,jj-1))^2+(r_b(2,jj)-r_b(2,jj-1))^2);
% 
% dH_b=side+front;
% end
% 
% if jj==length(r_a)
% side=-k(3)*(r_a(dim,jj)-r_b(dim,jj))*(sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2)-ll_s(jj))/sqrt((r_a(1,jj)-r_b(1,jj))^2+(r_a(2,jj)-r_b(2,jj))^2);
% %front=-k(2)*(r_b(dim,jj+1)-r_b(dim,jj))*(sqrt((r_b(1,jj+1)-r_b(1,jj))^2+(r_b(2,jj+1)-r_b(2,jj))^2)-ll_b(jj))/sqrt((r_b(1,jj+1)-r_b(1,jj))^2+(r_b(2,jj+1)-r_b(2,jj))^2);
% back=k(2)*(r_b(dim,jj)-r_b(dim,jj-1))*(sqrt((r_b(1,jj)-r_b(1,jj-1))^2+(r_b(2,jj)-r_b(2,jj-1))^2)-ll_b(jj-1))/sqrt((r_b(1,jj)-r_b(1,jj-1))^2+(r_b(2,jj)-r_b(2,jj-1))^2);
% 
% dH_b=side+back;               
% end
% bottom1=0;
%     bottom2=0;
% if dim==2 && jj==1
%     bottom1=gamma/lambda*exp(-1/2*(r_b(2,jj)+r_b(2,jj+1))/lambda);
% end  
% if dim==2 && jj==length(r_a)
%     bottom2=gamma/lambda*exp(-1/2*(r_b(2,jj)+r_b(2,jj-1))/lambda);
%     
%     %bottom1=-kappa/2*(rr(jj)-1/2*(r_b(1,jj)+r_b(1,jj+1)));
%     %bottom2=-kappa/2*(rr(jj-1)-1/2*(r_b(1,jj)+r_b(1,jj-1)));
% end
% 
% dH_b=dH_b+bottom1+bottom2;
% end

function dH_a=dH_a_func(ll_a,ll_s,r_a,r_b,kappa,jj,dim,k,g,delta_x,step,A_0,Area)
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

function dH_b=dH_b_func(ll_b,ll_s,r_a,r_b,kappa,jj,rr,dim,k,gamma,lambda, g, delta_x,A_0,Area)
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


bottom1=0;
bottom2=0;
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

function [r_a,r_b,rr,ll_a,ll_b,ll_s,ID,A_0,step,ID_initial,ID_edge] = ChainModelInitial(N, sigma,lmax,m) % Number of cells and noise in initial condition
    
    rng(0);
    r_a(1,:) = ((1+30:1:N+30)'-(N)/2);%  + rand(N,1)*sigma;% define cell center, initially at the 0.5, for example one cell would have center 0.5, then add noise
    
    rng(1)
    r_a(2,:)=ones(N,1);%+rand(N,1)*sigma;
 
    rng(2);
    r_b(1,:) = ((1+30:1:N+30)'-(N)/2);%  + rand(N,1)*sigma; % define cell center, initially at the 0.5, for example one cell would have center 0.5, then add noise
    
    
    rng(3)
    r_b(2,:)=zeros(N,1);%rand(N,1)*sigma
    
    % ID_initial marks the apical point of the boundary between N/E
    ID_initial=zeros(length(r_a),1);
    ID_initial((length(r_a)-1)/4+1)= 1;
    ID_initial((length(r_a)-1)*3/4+1)= 2;
    ID_initial(1)= -1;
    ID_initial(end)= -1;
    rr = 0.5*(r_b(1,2:end)+r_b(1,1:end-1));% + rand(N-1,1)*sigma; %attachment points (on the substrate)
    %l  = r(2:end)-r(1:end-1);%ones(N-1,1);
    %l = l-min(l)+.2;
    ll_a=1*ones(N-1,1); %ll_a = unifrnd(1,lmax,N-1,1)*4; %random numbers from distribution, to get size of cells %ecad top
    ll_a(length(ll_a)/4+1:length(ll_a)*3/4)=1*ones((N-1)/2,1);%ll_a(length(ll_a)/4+1:length(ll_a)*3/4)=unifrnd(1,lmax,(N-1)/2,1); %ncad top 0.1 0.5
    
    ll_b=1*ones(N-1,1);%ll_b = unifrnd(1,lmax,N-1,1); %random numbers from distribution, to get size of cells %ecad bottom
    ll_b(length(ll_b)/4+1:length(ll_b)*3/4)=1*ones((N-1)/2,1);%ll_b(length(ll_b)/4+1:length(ll_b)*3/4)=unifrnd(1,lmax,(N-1)/2,1); %ncad bottom
    
    ll_s=1*ones(N,1);%ll_b = unifrnd(1,lmax,N-1,1); %random numbers from distribution, to get size of cells %ecad bottom
    ll_s((length(ll_s)-1)/4+1:(length(ll_s)-1)*3/4+1)=1*ones((N-1)/2+1,1);
    ID=ones(length(ll_a),1);
    ID(length(ll_a)/4+1:length(ll_a)*3/4)= 2;
    
    A_0=ones(length(ll_a));
    
    step=0*ones(N-1,1);%ll_b = unifrnd(1,lmax,N-1,1); %random numbers from distribution, to get size of cells %ecad bottom
    step(length(ll_b)/4+1:length(ll_b)*3/4)=m*ones((N-1)/2,1);
    
    
    ID_edge=zeros(length(r_a),1);
    ID_edge((length(r_a)-1)/4+1)= 1;
    ID_edge((length(r_a)-1)*3/4+1)= 2;
   
end