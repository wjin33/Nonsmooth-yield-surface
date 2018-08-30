function [r_p,r_q,r_epsT,r_epsel,r_epsed,r_epsE,r_epsin,r_rho,r_f,r_Yd,r_energy] = DDM_CCP()
%%     Damage Consitutive Model (one element test?strain controlled test)
clc
close all

global nu0 E0 N_V0 a_0 Kc Ec Ko Eo alpha
global FTOL0 FTOL1 c0 c1 Iter

%% Orientation integral points and weights
Beta=33.2699078510*pi/180;
T1=cos(pi/4)*cos(pi/2-Beta);
T2=cos(Beta);
n_cosines=[1 0 0;
           0 1 0;
           0 0 1;
           sqrt(2)/2   sqrt(2)/2  0;
           sqrt(2)/2  -sqrt(2)/2  0;
           sqrt(2)/2   0          sqrt(2)/2;
           sqrt(2)/2   0         -sqrt(2)/2;
           0           sqrt(2)/2  sqrt(2)/2;
           0           sqrt(2)/2 -sqrt(2)/2;
           T1          T1         T2;
           T1          T1        -T2;
           T1         -T1         T2;
           T1         -T1        -T2;
           T1          T2         T1;
           T1          T2        -T1;
           T1         -T2         T1;
           T1         -T2        -T1;
           T2          T1         T1;
           T2          T1        -T1;
           T2         -T1         T1;
           T2         -T1        -T1;];
n_cosines(22:42,1:3)=-n_cosines(1:21,1:3);   
n_weight(1:3,1)=0.0265214244093;
n_weight(4:9,1)=0.0199301476312;
n_weight(10:21,1)=0.0250712367487;
n_weight(22:42,1)=n_weight(1:21,1);

%% Controlling Parameters for computation
FTOL0 = 1e-8; FTOL1=1e-20;
Iter=50;

%% Input and Store the material Parameters

% load('calied_DWCD8_13.mat');
% Pars(1,1:7) = Result'.*[10^10 1 10^2 10^-3 10^2 10^2 10^-5];
% load('calied_DWCD8_tension.mat');
% Pars(1,8:9)=Result';

% Pars =[3.85E10 0.34 46 0.002 36 32 1.65E-5 18 23];

Pars = [5.35E10 0.35 960 0.05 278.9 116.6 10^-5 35.9 20.6];
E0   = Pars(1);
nu0  = Pars(2);
N_V0 = Pars(3);
a_0  = Pars(4);
Kc   = Pars(5);
Ec   = Pars(6);
alpha= Pars(7);
Ko   = Pars(8);
Eo   = Pars(9);

c0=(16/3)*(1-nu0^2)/E0;
c1=(32/3)*(1-nu0^2)/(2-nu0)/E0;


%% Load path illustration:
%  Eta = %%  1 = iso, 2 = uniaxial, 3 = shear
% steps -- how many loading steps for this simulation, normally its the size of Eta
% chargEps1 -- applied strain component for each increments(Only strain controlled tests are coded as stress controll do not nead to iteration)
% ninc -- Increments used for each loading step


%% Load path example
%uniaxial-tensile-compressive
% Eta = [2 2 2];
% steps = 3;
% chargEps1 = [-0.001 0.0005 0.02];
% ninc = [100 10 200];


% % cyclic-uniaxial-compressive
% Eta = [2 2 2 2];
% steps = 4;
% chargEps1 = [0.02 -0.005 0.01 -0.005];
% ninc = [200 10 200 1];


%pure shear 
Eta = [3];
steps = 1;
chargEps1 = [0.01];
ninc = [200];

%traxial compression
% Eta = [1 2];
% steps = 2;
% chargEps1 = [0.001 0.02];
% ninc = [20 200];



%%  INITIALIZATIONS:
sigmaT = zeros(3,3); % non-loading
sigmaT_v = mat2_mat1(sigmaT);
%========================================================================
%         Strain decomposition:
%         epsT = epsel + epsed + epsin
%========================================================================
epsT = zeros(3,3);
epsel = zeros(3,3);
epsed = zeros(3,3);
epsE = zeros(3,3);
epsin = zeros(3,3);
Yd = zeros(42,1);
rho_0= N_V0*a_0^3*ones(42,1);
rho=rho_0;
energy = zeros(5,1);
[~,matS0] = matDO1(Pars,rho_0,n_cosines,n_weight,sigmaT);
%% transfer matrix to vector
epsT_v = mat2_mat1(epsT);
epsel_v = mat2_mat1(epsel);
epsed_v = mat2_mat1(epsed);
epsE_v = mat2_mat1(epsE);
epsin_v = mat2_mat1(epsin);
%% Storage
r_p(1,:) = sum(sigmaT_v)/3;
r_q(1,:) = sigmaT(1,1)-sigmaT(2,2);
r_sigmaT(1,:) = sigmaT_v;
r_epsel(1,:) = epsel_v;
r_epsed(1,:) = epsed_v;
r_epsin(1,:) = epsin_v;
r_epsE(1,:) = epsE_v;
r_epsT(1,:) = epsT_v;
r_Yd(1,:) = Yd';
r_rho(1,:)=rho';
r_energy(1,:) = energy';
% yield function
[r_f(1,:),r_Yd(1,:),~] = fdDP(Pars,rho,n_cosines,n_weight,sigmaT); %trial test

%% SIMULATION (for a single stress path component):


% tinc = 2;
for icharg = 1:steps 
    disp(['============= load step #',num2str(icharg),' ============='])
    
    for inc = 1:ninc(icharg) % load increments
        disp(['              increments #',num2str(inc),'              '])
                
        [matDOm,~] = matDO1(Pars,rho,n_cosines,n_weight,sigmaT);% matDOm-stiffness matrix 4th order

        if Eta(icharg) == 1 %iso
            deps = chargEps1(icharg)/ninc(icharg)*[1 0 0;0 1 0;0 0 1;];
            dsig = Aijkl_Bkl(matDOm,deps);
        elseif Eta(icharg) == 2 %uniaxial 
            deps = chargEps1(icharg)/ninc(icharg)*[1 0 0;0 0 0;0 0 0;];
            dsig = Aijkl_Bkl(matDOm,deps);
        elseif Eta(icharg) == 3 % shear
            deps = chargEps1(icharg)/ninc(icharg)*[0 1 0;1 0 0;0 0 0];
            dsig = Aijkl_Bkl(matDOm,deps);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        sigmaTri = sigmaT + dsig;
        [fd,Yd,activated] = fdDP(Pars,rho,n_cosines,n_weight,sigmaTri); %trial test                 
% fd
% activated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        if fd <= FTOL0 % elastic increment
            
            % Elastic trial:
            depsel = Aijkl_Bkl(matS0,dsig);
            depsin = zeros(3,3);
            depsE = deps;
            depsed=depsE-depsel;
            drho=zeros(42,1);
            depsT= depsE + depsin;

            sigmaT = sigmaT + dsig;
            sigEner = sigmaT-0.5*dsig;
            epsE = epsE + depsE;
            epsT = epsT + depsT ;
            epsin = epsin + depsin;
            epsel = epsel + depsel;
            epsed = epsed + depsed;
            rho = rho + drho;

            pn = trace(sigmaT)/3;
            qn = sigmaT(1,1) - sigmaT(2,2);
            % transfer matrix to vector
            sigmaT_v = mat2_mat1(sigmaT);
            epsT_v = mat2_mat1(epsT);
            epsE_v = mat2_mat1(epsE);
            epsel_v = mat2_mat1(epsel); 
            epsed_v = mat2_mat1(epsed);
            epsin_v = mat2_mat1(epsin);

            energy(1,1)=Aij_Bij(sigEner,depsT);
            energy(2,1)=Aij_Bij(sigEner,depsE);
            energy(3,1)=Aij_Bij(sigEner,depsel);
            energy(4,1)=Aij_Bij(sigEner,depsin);
            energy(5,1)= 0 ; %No damage accumulation
            
            %Start to store components
            end_r = length(r_p)+1;
            r_p(end_r) = pn;
            r_q(end_r) = qn;
            r_sigmaT(end_r,:)=sigmaT_v;
            r_epsT(end_r,:) = epsT_v;
            r_epsel(end_r,:) = epsel_v;
            r_epsed(end_r,:) = epsed_v;
            r_epsE(end_r,:) = epsE_v;
            r_epsin(end_r,:) = epsin_v;
            r_Yd(end_r,:) = Yd';
            r_rho(end_r,:) = rho';
            r_f(end_r,:) = fd';
            r_energy(end_r,:) =r_energy(end_r-1,1:5)+energy';
            
        else % inelastic increment
            
            incinc=1;
            flag=0;
            while flag ~= 1
                [~,n]=size(activated);
                sigmaT_i=sigmaTri;
                drho_i=zeros(n,1);
                R_i=zeros(3,3);
                depsin_i=zeros(3,3);
                while incinc<Iter
                    [ddepsin,ddrho,dsig] = fd_lam(Pars, rho, drho_i, sigmaT_i,n_cosines,n_weight,activated,R_i);
                    depsin_f=depsin_i+ddepsin;
                    sigmaT_f=sigmaT_i+dsig;
                    drho_f=drho_i+ddrho;                 
                    detec=0;
                    for i=1:n
                        if drho_f(i,1)<0
                            activated(1,i:n-1)=activated(1,i+1:n);
                            detec=detec+1;
                        end
                    end
                    if detec ~= 0
                        activated_new=activated(1,1:n-detec);
                        activated=activated_new;
                        break
                    end
                    
                    depsin_check=zeros(3,3);
                    for i=1:n
                        i_point=activated(1,i);
                        dg_dsig =  Dg_Dsigma(Pars, n_cosines(i_point,:), n_weight(i_point,1), sigmaT_f);
                        depsin_check = depsin_check + drho_f(i)*dg_dsig;
                    end

                    R_f = -depsin_f+depsin_check;
                    Res = abs(trace(R_f))/3;
                    
                    [fd_all] = fdDP1(Pars,rho,drho_f,n_cosines,n_weight,sigmaT_f,activated);
                    
                    if (fd_all < FTOL0) && (Res < FTOL1)
                        flag=1;
                        break % Iteration satisfy requirements, Done
                    else
                        drho_i=drho_f;
                        sigmaT_i=sigmaT_f;
                        depsin_i=depsin_f;
                        R_i=R_f;
                        incinc=incinc+1;
                    end
                end
            end
            for i=1:n
                i_point=activated(1,i);
                rho(i_point,1)=rho(i_point,1)+drho_f(i,1);
            end
            
            dsig=sigmaT_f-sigmaT;
            depsin=depsin_f;
            depsT= deps;
            depsel = Aijkl_Bkl(matS0,dsig);
            depsed = deps  - depsel-depsin;
            depsE = depsel + depsed;
            
            
            sigmaT = sigmaT_f;
            sigEner = sigmaT-0.5*dsig;
            epsE = epsE + depsE;
            epsT = epsT + depsT ;
            epsin = epsin + depsin;
            epsel = epsel + depsel;
            epsed = epsed + depsed;

            pn = trace(sigmaT)/3;
            qn = sigmaT(1,1) - sigmaT(2,2);
            % transfer matrix to vector
            sigmaT_v = mat2_mat1(sigmaT);
            epsT_v = mat2_mat1(epsT);
            epsE_v = mat2_mat1(epsE);
            epsel_v = mat2_mat1(epsel); 
            epsed_v = mat2_mat1(epsed);
            epsin_v = mat2_mat1(epsin);
            energy(1,1)=Aij_Bij(sigEner,depsT);
            energy(2,1)=Aij_Bij(sigEner,depsE);
            energy(3,1)=Aij_Bij(sigEner,depsel);
            energy(4,1)=Aij_Bij(sigEner,depsin);
            [fd,Yd,~] = fdDP(Pars,rho,n_cosines,n_weight,sigmaT);
            for i=1:n
                i_point=activated(1,i);
                energy(5,1) = energy(5,1) + Yd(i_point,1)*drho_f(i,1);
            end
            
            %Start to store components
            end_r = length(r_p)+1;
            r_p(end_r) = pn;
            r_q(end_r) = qn;
            r_sigmaT(end_r,:)=sigmaT_v;
            r_epsT(end_r,:) = epsT_v;
            r_epsel(end_r,:) = epsel_v;
            r_epsed(end_r,:) = epsed_v;
            r_epsE(end_r,:) = epsE_v;
            r_epsin(end_r,:) = epsin_v;
            r_Yd(end_r,:) = Yd';
            r_rho(end_r,:) = rho';
            r_f(end_r,:) = fd';
            r_energy(end_r,:) =r_energy(end_r-1,1:5)+energy';                                           
        end
    end % increments

end % loading parts


% save confined_compression_200  r_sigmaT r_epsT r_rho
% save shear_200  r_sigmaT r_epsT r_rho

% return

figure('Name','q(eps1^T,eps3^T)','NumberTitle','off');
plot(r_epsT(:,1),r_q/1000000,'-b',r_epsT(:,3),r_q/1000000,'-r','Linewidth',2)
xlabel('Strain, \epsilon^T','FontSize',20)
ylabel('Deviatoric stress, \sigma_1-\sigma_3 (MPa)','FontSize',20)
legend('\epsilon_1^T', '\epsilon_3^T','Location','Best')
grid
set(gca,'FontSize',20)
set(gcf,'Units','inches')
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);    
% print(gcf,'-dpdf','DDM_kappa_1-3-strain.pdf')



figure('Name','\tau-q','NumberTitle','off');
plot(r_epsT(:,4),r_sigmaT(:,4)/1000000,'-b','Linewidth',2)
xlabel('Strain, \epsilon^T_{12}','FontSize',20)
ylabel('Deviatoric stress, \tau_{12} (MPa)','FontSize',20)
% legend('\epsilon_1^T', '\epsilon_3^T','Location','Best')
grid
set(gca,'FontSize',20)
set(gcf,'Units','inches')
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);    
% print(gcf,'-dpdf','DWCD8_gamma12_sigma12-strain.pdf')






% figure('Name','q(eps1^{ed},eps3^{ed})','NumberTitle','off');
% plot(r_epsed(:,1),r_q/1000000,'-b',r_epsed(:,3),r_q/1000000,'-r','Linewidth',2)
% xlabel('Axial strain, \epsilon^{ed}','FontSize',20)
% ylabel('Deviatoric stress, q (MPa)','FontSize',20)
% legend('\epsilon_1^{ed}', '\epsilon_3^{ed}','Location','Best')
% grid
% set(gca,'FontSize',20)
% print -depsc 'DWCD_q_eps_ed_1-3.eps'
% 
% 
% figure('Name','q(eps1^{id},eps3^{id})','NumberTitle','off');
% plot(r_epsin(:,1),r_q/1000000,'-b',r_epsin(:,3),r_q/1000000,'-r','Linewidth',2)
% xlabel('Axial strain, \epsilon^{id}','FontSize',20)
% ylabel('Deviatoric stress, q (MPa)','FontSize',20)
% legend('\epsilon_1^{id}', '\epsilon_3^{id}','Location','Best')
% grid
% set(gca,'FontSize',20)
% print -depsc 'DWCD_q_eps_id_1-3.eps'



figure
plot3(r_rho(:,1)*n_cosines(1,2), r_rho(:,1)*n_cosines(1,3), r_rho(:,1)*n_cosines(1,1),'Linewidth',2)
hold on
plot3(r_rho(:,2)*n_cosines(2,2), r_rho(:,2)*n_cosines(2,3), r_rho(:,2)*n_cosines(2,1),'Linewidth',2)
plot3(r_rho(:,3)*n_cosines(3,2), r_rho(:,3)*n_cosines(3,3), r_rho(:,3)*n_cosines(3,1),'Linewidth',2)
plot3(r_rho(:,4)*n_cosines(4,2), r_rho(:,4)*n_cosines(4,3), r_rho(:,4)*n_cosines(4,1),'Linewidth',2)
plot3(r_rho(:,5)*n_cosines(5,2), r_rho(:,5)*n_cosines(5,3), r_rho(:,5)*n_cosines(5,1),'Linewidth',2)
plot3(r_rho(:,6)*n_cosines(6,2), r_rho(:,6)*n_cosines(6,3), r_rho(:,6)*n_cosines(6,1),'Linewidth',2)
plot3(r_rho(:,7)*n_cosines(7,2), r_rho(:,7)*n_cosines(7,3), r_rho(:,7)*n_cosines(7,1),'Linewidth',2)
plot3(r_rho(:,8)*n_cosines(8,2), r_rho(:,8)*n_cosines(8,3), r_rho(:,8)*n_cosines(8,1),'Linewidth',2)
plot3(r_rho(:,9)*n_cosines(9,2), r_rho(:,9)*n_cosines(9,3), r_rho(:,9)*n_cosines(9,1),'Linewidth',2)
plot3(r_rho(:,10)*n_cosines(10,2), r_rho(:,10)*n_cosines(10,3), r_rho(:,10)*n_cosines(10,1),'Linewidth',2)
plot3(r_rho(:,11)*n_cosines(11,2), r_rho(:,11)*n_cosines(11,3), r_rho(:,11)*n_cosines(11,1),'Linewidth',2)
plot3(r_rho(:,12)*n_cosines(12,2), r_rho(:,12)*n_cosines(12,3), r_rho(:,12)*n_cosines(12,1),'Linewidth',2)
plot3(r_rho(:,13)*n_cosines(13,2), r_rho(:,13)*n_cosines(13,3), r_rho(:,13)*n_cosines(13,1),'Linewidth',2)
plot3(r_rho(:,14)*n_cosines(14,2), r_rho(:,14)*n_cosines(14,3), r_rho(:,14)*n_cosines(14,1),'Linewidth',2)
plot3(r_rho(:,15)*n_cosines(15,2), r_rho(:,15)*n_cosines(15,3), r_rho(:,15)*n_cosines(15,1),'Linewidth',2)
plot3(r_rho(:,16)*n_cosines(16,2), r_rho(:,16)*n_cosines(16,3), r_rho(:,16)*n_cosines(16,1),'Linewidth',2)
plot3(r_rho(:,17)*n_cosines(17,2), r_rho(:,17)*n_cosines(17,3), r_rho(:,17)*n_cosines(17,1),'Linewidth',2)
plot3(r_rho(:,18)*n_cosines(18,2), r_rho(:,18)*n_cosines(18,3), r_rho(:,18)*n_cosines(18,1),'Linewidth',2)
plot3(r_rho(:,19)*n_cosines(19,2), r_rho(:,19)*n_cosines(19,3), r_rho(:,19)*n_cosines(19,1),'Linewidth',2)
plot3(r_rho(:,20)*n_cosines(20,2), r_rho(:,20)*n_cosines(20,3), r_rho(:,20)*n_cosines(20,1),'Linewidth',2)
plot3(r_rho(:,21)*n_cosines(21,2), r_rho(:,21)*n_cosines(21,3), r_rho(:,21)*n_cosines(21,1),'Linewidth',2)
plot3(r_rho(:,22)*n_cosines(22,2), r_rho(:,22)*n_cosines(22,3), r_rho(:,22)*n_cosines(22,1),'Linewidth',2)
plot3(r_rho(:,23)*n_cosines(23,2), r_rho(:,23)*n_cosines(23,3), r_rho(:,23)*n_cosines(23,1),'Linewidth',2)
plot3(r_rho(:,24)*n_cosines(24,2), r_rho(:,24)*n_cosines(24,3), r_rho(:,24)*n_cosines(24,1),'Linewidth',2)
plot3(r_rho(:,25)*n_cosines(25,2), r_rho(:,25)*n_cosines(25,3), r_rho(:,25)*n_cosines(25,1),'Linewidth',2)
plot3(r_rho(:,26)*n_cosines(26,2), r_rho(:,26)*n_cosines(26,3), r_rho(:,26)*n_cosines(26,1),'Linewidth',2)
plot3(r_rho(:,27)*n_cosines(27,2), r_rho(:,27)*n_cosines(27,3), r_rho(:,27)*n_cosines(27,1),'Linewidth',2)
plot3(r_rho(:,28)*n_cosines(28,2), r_rho(:,28)*n_cosines(28,3), r_rho(:,28)*n_cosines(28,1),'Linewidth',2)
plot3(r_rho(:,29)*n_cosines(29,2), r_rho(:,29)*n_cosines(29,3), r_rho(:,29)*n_cosines(29,1),'Linewidth',2)
plot3(r_rho(:,30)*n_cosines(30,2), r_rho(:,30)*n_cosines(30,3), r_rho(:,30)*n_cosines(30,1),'Linewidth',2)
plot3(r_rho(:,31)*n_cosines(31,2), r_rho(:,31)*n_cosines(31,3), r_rho(:,31)*n_cosines(31,1),'Linewidth',2)
plot3(r_rho(:,32)*n_cosines(32,2), r_rho(:,32)*n_cosines(32,3), r_rho(:,32)*n_cosines(32,1),'Linewidth',2)
plot3(r_rho(:,33)*n_cosines(33,2), r_rho(:,33)*n_cosines(33,3), r_rho(:,33)*n_cosines(33,1),'Linewidth',2)
plot3(r_rho(:,34)*n_cosines(34,2), r_rho(:,34)*n_cosines(34,3), r_rho(:,34)*n_cosines(34,1),'Linewidth',2)
plot3(r_rho(:,35)*n_cosines(35,2), r_rho(:,35)*n_cosines(35,3), r_rho(:,35)*n_cosines(35,1),'Linewidth',2)
plot3(r_rho(:,36)*n_cosines(36,2), r_rho(:,36)*n_cosines(36,3), r_rho(:,36)*n_cosines(36,1),'Linewidth',2)
plot3(r_rho(:,37)*n_cosines(37,2), r_rho(:,37)*n_cosines(37,3), r_rho(:,37)*n_cosines(37,1),'Linewidth',2)
plot3(r_rho(:,38)*n_cosines(38,2), r_rho(:,38)*n_cosines(38,3), r_rho(:,38)*n_cosines(38,1),'Linewidth',2)
plot3(r_rho(:,39)*n_cosines(39,2), r_rho(:,39)*n_cosines(39,3), r_rho(:,39)*n_cosines(39,1),'Linewidth',2)
plot3(r_rho(:,40)*n_cosines(40,2), r_rho(:,40)*n_cosines(40,3), r_rho(:,40)*n_cosines(40,1),'Linewidth',2)
plot3(r_rho(:,41)*n_cosines(41,2), r_rho(:,41)*n_cosines(41,3), r_rho(:,41)*n_cosines(41,1),'Linewidth',2)
plot3(r_rho(:,42)*n_cosines(42,2), r_rho(:,42)*n_cosines(42,3), r_rho(:,42)*n_cosines(42,1),'Linewidth',2)
xlabel('\rho_2 ','FontSize',15)
ylabel('\rho_3 ','FontSize',15)
zlabel('Damage density, \rho_1 ','FontSize',15)
grid on
axis equal
set(gca,'FontSize',15)
set(gcf,'Units','inches')
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);    
% print(gcf,'-dpdf','DDM_3D_density-odometer.pdf')



figure
plot(r_q/1000000,r_rho(:,1),'b-','Linewidth',3)
hold on
plot(r_q/1000000,r_rho(:,4),'r--','Linewidth',3)
plot(r_q/1000000,r_rho(:,10),'k-','Linewidth',3)
plot(r_q/1000000,r_rho(:,18),'g-','Linewidth',3)
legend('\rho_1','\rho_4','\rho_{10}','\rho_{18}','Location','Northwest')
xlabel('Axial stess, \sigma_1 (MPa)','FontSize',20)
ylabel('Crack density, \rho','FontSize',20)
grid
set(gca,'FontSize',20)
set(gcf,'Units','inches')
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);    
% print(gcf,'-dpdf','DWCD8_2D_q_density-strain.pdf')



figure
plot(r_epsT(:,4),r_rho(:,1),'b-','Linewidth',3)
hold on
plot(r_epsT(:,4),r_rho(:,4),'r--','Linewidth',3)
plot(r_epsT(:,4),r_rho(:,10),'k-','Linewidth',3)
plot(r_epsT(:,4),r_rho(:,18),'g-','Linewidth',3)
legend('\rho_1','\rho_4','\rho_{10}','\rho_{18}','Location','Northwest')
xlabel('Shear strain, \epsilon_{12}','FontSize',20)
ylabel('Crack density, \rho','FontSize',20)
grid
set(gca,'FontSize',20)
set(gcf,'Units','inches')
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);    
% print(gcf,'-dpdf','DWCD8_2D_gamma12_density-strain.pdf')



figure1 = figure('NumberTitle','off','Name','Energy-tension');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot4 = plot(r_sigmaT(:,1)/1000000,r_energy(:,1),'Parent',axes1,'LineWidth',2);
plot5 = plot(r_sigmaT(:,1)/1000000,r_energy(:,2),'Parent',axes1,'Linewidth',2);
plot6 = plot(r_sigmaT(:,1)/1000000,r_energy(:,4),'Parent',axes1,'Linewidth',2);
% plot7 = plot(r_sigmaT(:,1)/1000000,r_energy(:,5),'Parent',axes1,'Linewidth',3);

set(plot4,'LineStyle','-','Color',[0 0 1],'DisplayName','External work');
set(plot5,'LineStyle','--','Color',[1 0 0],'DisplayName','Elastic strain energy');
set(plot6,'LineStyle','-.','Color',[0 0 0.5],'DisplayName','Inelastic strain energy');
% set(plot7,'LineStyle','-.','Color',[0 1 1],'DisplayName','Crack debonding');
% Create xlabel
xlabel('Axial stress, \sigma_{1} (MPa)','FontSize',20);
% Create ylabel
ylabel('Energy (J)','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','Best','FontSize',15);
set(gcf,'Units','inches')
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);    
% print(gcf,'-dpdf','DWCD8_energy-strain.pdf')

end


%==========================================================================
%
% ----------------------------SUB-FUNCTIOIN--------------------------------
%
%==========================================================================


function [ scalar ] = Aij_Bij( A,B )
%   Double contraction: scalar = A_ij*B_ij

scalar = 0;
for i= 1:3
    for j=1:3 
        scalar = scalar + A(i,j)*B(i,j); 
    end
end
   
end

function [ C ] = Aijkl_Bij( A,B )

C = zeros(3,3);
for k = 1:3
    for l = 1:3
        for i = 1:3
            for j = 1:3
              C(k,l) = C(k,l)+A(i,j,k,l)*B(i,j);
            end
        end
    end
end

end

function [ C ] = Aijkl_Bkl( A,B )

C = zeros(3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
              C(i,j) = C(i,j)+A(i,j,k,l)*B(k,l);
            end
        end
    end
end

end

function [ vector ] = mat2_mat1( matrix )
%==========================================================================
%
%    MAT2_MAT1 
%    Transfer a 3*3 matrix to 6*1 vextor (default format in ABAQUS)
%
%    sig_11 sig_12 sig_13
%    sig_21 sig_22 sig_23 ==> [sig_11 sig_22 sig_33 sig_12 sig_13 sig_23]^T
%    sig_31 sig_32 sig_33
%
%==========================================================================
vector = zeros(1,6);
for i = 1:3
    vector(i) = matrix(i,i);
    for j = i+1:3
        vector(i+j+1) = matrix(i,j);
    end
end

end

function [ tensor ] = mat2_mat4( matrix,coe )
      tensor(1:3,1:3,1:3,1:3) = 0;
      if coe==1
          coe1=1;
          coe2=1;
      elseif coe==2
          coe1=2;
          coe2=4;
      end
      tensor(1,1,1,1) = matrix(1,1);
      tensor(1,1,2,2) = matrix(1,2);
      tensor(1,1,3,3) = matrix(1,3);
      tensor(1,1,1,2) = matrix(1,4);
      tensor(1,1,2,1) = matrix(1,4)/coe1;
      tensor(1,1,2,3) = matrix(1,5)/coe1;
      tensor(1,1,3,2) = matrix(1,5)/coe1;
      tensor(1,1,1,3) = matrix(1,6)/coe1;
      tensor(1,1,3,1) = matrix(1,6)/coe1;

      tensor(2,2,1,1) = matrix(2,1);
      tensor(2,2,2,2) = matrix(2,2);
      tensor(2,2,3,3) = matrix(2,3);
      tensor(2,2,1,2) = matrix(2,4)/coe1;
      tensor(2,2,2,1) = matrix(2,4)/coe1;
      tensor(2,2,2,3) = matrix(2,5)/coe1;
      tensor(2,2,3,2) = matrix(2,5)/coe1;
      tensor(2,2,1,3) = matrix(2,6)/coe1;
      tensor(2,2,3,1) = matrix(2,6)/coe1;

      tensor(3,3,1,1) = matrix(3,1);
      tensor(3,3,2,2) = matrix(3,2);
      tensor(3,3,3,3) = matrix(3,3);
      tensor(3,3,1,2) = matrix(3,4)/coe1;
      tensor(3,3,2,1) = matrix(3,4)/coe1;
      tensor(3,3,2,3) = matrix(3,5)/coe1;
      tensor(3,3,3,2) = matrix(3,5)/coe1;
      tensor(3,3,1,3) = matrix(3,6)/coe1;
      tensor(3,3,3,1) = matrix(3,6)/coe1;

      tensor(1,2,1,1) = matrix(4,1)/coe1;
      tensor(1,2,2,2) = matrix(4,2)/coe1;
      tensor(1,2,3,3) = matrix(4,3)/coe1;
      tensor(1,2,1,2) = matrix(4,4)/coe2;
      tensor(1,2,2,1) = matrix(4,4)/coe2;
      tensor(1,2,2,3) = matrix(4,5)/coe2;
      tensor(1,2,3,2) = matrix(4,5)/coe2;
      tensor(1,2,1,3) = matrix(4,6)/coe2;
      tensor(1,2,3,1) = matrix(4,6)/coe2;

      tensor(2,3,1,1) = matrix(5,1)/coe1;
      tensor(2,3,2,2) = matrix(5,2)/coe1;
      tensor(2,3,3,3) = matrix(5,3)/coe1;
      tensor(2,3,1,2) = matrix(5,4)/coe2;
      tensor(2,3,2,1) = matrix(5,4)/coe2;
      tensor(2,3,2,3) = matrix(5,5)/coe2;
      tensor(2,3,3,2) = matrix(5,5)/coe2;
      tensor(2,3,1,3) = matrix(5,6)/coe2;
      tensor(2,3,3,1) = matrix(5,6)/coe2;

      tensor(1,3,1,1) = matrix(6,1)/coe1;
      tensor(1,3,2,2) = matrix(6,2)/coe1;
      tensor(1,3,3,3) = matrix(6,3)/coe1;
      tensor(1,3,1,2) = matrix(6,4)/coe2;
      tensor(1,3,2,1) = matrix(6,4)/coe2;
      tensor(1,3,2,3) = matrix(6,5)/coe2;
      tensor(1,3,3,2) = matrix(6,5)/coe2;
      tensor(1,3,1,3) = matrix(6,6)/coe2;
      tensor(1,3,3,1) = matrix(6,6)/coe2;
      
      tensor(2,1,1,1) = matrix(4,1)/coe1;
      tensor(2,1,2,2) = matrix(4,2)/coe1;
      tensor(2,1,3,3) = matrix(4,3)/coe1;
      tensor(2,1,1,2) = matrix(4,4)/coe2;
      tensor(2,1,2,1) = matrix(4,4)/coe2;
      tensor(2,1,2,3) = matrix(4,5)/coe2;
      tensor(2,1,3,2) = matrix(4,5)/coe2;
      tensor(2,1,1,3) = matrix(4,6)/coe2;
      tensor(2,1,3,1) = matrix(4,6)/coe2;

      tensor(3,2,1,1) = matrix(5,1)/coe1;
      tensor(3,2,2,2) = matrix(5,2)/coe1;
      tensor(3,2,3,3) = matrix(5,3)/coe1;
      tensor(3,2,1,2) = matrix(5,4)/coe2;
      tensor(3,2,2,1) = matrix(5,4)/coe2;
      tensor(3,2,2,3) = matrix(5,5)/coe2;
      tensor(3,2,3,2) = matrix(5,5)/coe2;
      tensor(3,2,1,3) = matrix(5,6)/coe2;
      tensor(3,2,3,1) = matrix(5,6)/coe2;

      tensor(3,1,1,1) = matrix(6,1)/coe1;
      tensor(3,1,2,2) = matrix(6,2)/coe1;
      tensor(3,1,3,3) = matrix(6,3)/coe1;
      tensor(3,1,1,2) = matrix(6,4)/coe2;
      tensor(3,1,2,1) = matrix(6,4)/coe2;
      tensor(3,1,2,3) = matrix(6,5)/coe2;
      tensor(3,1,3,2) = matrix(6,5)/coe2;
      tensor(3,1,1,3) = matrix(6,6)/coe2;
      tensor(3,1,3,1) = matrix(6,6)/coe2;

end

function [ matrix ] = mat4_mat2( tensor,coe )

      matrix = zeros (6,6);
      if coe==1
          coe1=1;
          coe2=1;
      elseif coe==2
          coe1=2;
          coe2=4;
      end
      matrix(1,1)=tensor(1,1,1,1);
      matrix(1,2)=tensor(1,1,2,2);
      matrix(1,3)=tensor(1,1,3,3);
      matrix(1,4)=tensor(1,1,1,2)*coe1;
      matrix(1,5)=tensor(1,1,2,3)*coe1;
      matrix(1,6)=tensor(1,1,1,3)*coe1;

      matrix(2,1)=tensor(2,2,1,1);
      matrix(2,2)=tensor(2,2,2,2);
      matrix(2,3)=tensor(2,2,3,3);
      matrix(2,4)=tensor(2,2,1,2)*coe1;
      matrix(2,5)=tensor(2,2,2,3)*coe1;
      matrix(2,6)=tensor(2,2,1,3)*coe;

      matrix(3,1)=tensor(3,3,1,1);
      matrix(3,2)=tensor(3,3,2,2);
      matrix(3,3)=tensor(3,3,3,3);
      matrix(3,4)=tensor(3,3,1,2)*coe1;
      matrix(3,5)=tensor(3,3,2,3)*coe1;
      matrix(3,6)=tensor(3,3,1,3)*coe1;

      matrix(4,1)=tensor(1,2,1,1)*coe1;
      matrix(4,2)=tensor(1,2,2,2)*coe1;
      matrix(4,3)=tensor(1,2,3,3)*coe1;
      matrix(4,4)=tensor(1,2,1,2)*coe2;
      matrix(4,5)=tensor(1,2,2,3)*coe2;
      matrix(4,6)=tensor(1,2,1,3)*coe2;

      matrix(5,1)=tensor(2,3,1,1)*coe1;
      matrix(5,2)=tensor(2,3,2,2)*coe1;
      matrix(5,3)=tensor(2,3,3,3)*coe1;
      matrix(5,4)=tensor(2,3,1,2)*coe2;
      matrix(5,5)=tensor(2,3,2,3)*coe2;
      matrix(5,6)=tensor(2,3,1,3)*coe2;

      matrix(6,1)=tensor(1,3,1,1)*coe1;
      matrix(6,2)=tensor(1,3,2,2)*coe1;
      matrix(6,3)=tensor(1,3,3,3)*coe1;
      matrix(6,4)=tensor(1,3,1,2)*coe2;
      matrix(6,5)=tensor(1,3,2,3)*coe2;
      matrix(6,6)=tensor(1,3,1,3)*coe2;
end

function [matDz, matS] = matDO1(Pars,rho,n_cosines,n_weight,sigmaT)

% global E0 nu0 c0 c1

E0=Pars(1);
nu0=Pars(2);
c0=(16/3)*(1-nu0^2)/E0;
c1=(32/3)*(1-nu0^2)/(2-nu0)/E0;
b1 = (1+nu0)/E0/2;
b2 = nu0/E0;
E = eye(3);

NN(1:3,1:3,1:3,1:3) = 0;
Tri(1:3,1:3,1:3,1:3) = 0;
for i_point=1:42
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    if n_cosines(i_point,:)*sigmaT*n_cosines(i_point,:)'<=0
                        NN(i,j,k,l) = NN(i,j,k,l) + rho(i_point,1)*n_weight(i_point,1)*c0*n_cosines(i_point,i)*n_cosines(i_point,j)...
                                  *n_cosines(i_point,k)*n_cosines(i_point,l);
                    end
                    Tri(i,j,k,l) = Tri(i,j,k,l) + rho(i_point,1)*n_weight(i_point,1)*c1*(0.25*(n_cosines(i_point,i)*n_cosines(i_point,k)*E(j,l)+...
                                   n_cosines(i_point,i)*n_cosines(i_point,l)*E(j,k)+E(i,k)*n_cosines(i_point,j)*n_cosines(i_point,l)...
                                   +E(i,l)*n_cosines(i_point,j)*n_cosines(i_point,k))-n_cosines(i_point,i)*n_cosines(i_point,j)...
                                  *n_cosines(i_point,k)*n_cosines(i_point,l));
                end
            end
        end
    end
end
    
matS0(1:3,1:3,1:3,1:3) = 0;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                matS0(i,j,k,l) = b1*(E(i,k)*E(j,l)+E(i,l)*E(j,k))-b2*E(i,j)*E(k,l);
            end
        end
    end
end    
    
matS=matS0+NN+Tri;
   
matS_2 = mat4_mat2(matS,2);
E6 = eye(6);
matDz_2 = matS_2\E6;
matDz = mat2_mat4(matDz_2,1);

end

function [matDz, matS] = matDO2(Pars,rho,drho,activated,n_cosines,n_weight,sigmaT)

% global E0 nu0 c0 c1

E0=Pars(1);
nu0=Pars(2);
c0=(16/3)*(1-nu0^2)/E0;
c1=(32/3)*(1-nu0^2)/(2-nu0)/E0;
b1 = (1+nu0)/E0/2;
b2 = nu0/E0;
E = eye(3);

[~,n]=size(activated);
for i=1:n
    i_point=activated(1,i);
    rho(i_point,1)=rho(i_point,1)+drho(i,1);
end

NN(1:3,1:3,1:3,1:3) = 0;
Tri(1:3,1:3,1:3,1:3) = 0;
for i_point=1:42
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    if n_cosines(i_point,:)*sigmaT*n_cosines(i_point,:)'<=0
                        NN(i,j,k,l) = NN(i,j,k,l) + rho(i_point,1)*n_weight(i_point,1)*c0*n_cosines(i_point,i)*n_cosines(i_point,j)...
                                  *n_cosines(i_point,k)*n_cosines(i_point,l);
                    end
                    Tri(i,j,k,l) = Tri(i,j,k,l) + rho(i_point,1)*n_weight(i_point,1)*c1*(0.25*(n_cosines(i_point,i)*n_cosines(i_point,k)*E(j,l)+...
                                   n_cosines(i_point,i)*n_cosines(i_point,l)*E(j,k)+E(i,k)*n_cosines(i_point,j)*n_cosines(i_point,l)...
                                   +E(i,l)*n_cosines(i_point,j)*n_cosines(i_point,k))-n_cosines(i_point,i)*n_cosines(i_point,j)...
                                  *n_cosines(i_point,k)*n_cosines(i_point,l));
                end
            end
        end
    end
end
    
matS0(1:3,1:3,1:3,1:3) = 0;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                matS0(i,j,k,l) = b1*(E(i,k)*E(j,l)+E(i,l)*E(j,k))-b2*E(i,j)*E(k,l);
            end
        end
    end
end    


matS=matS0+NN+Tri;
   
matS_2 = mat4_mat2(matS,2);
E6 = eye(6);
matDz_2 = matS_2\E6;
matDz = mat2_mat4(matDz_2,1);

end

function [fd,Yd,activated] = fdDP(Pars,rho,n_cosines,n_weight,sigmaT)
%========================================================
% determine yield criteria with activated carck families
%========================================================

activated=[];



fd=zeros(42,1);
Yd=zeros(42,1);

i=1;
for i_point=1:42
    [Yd(i_point,1),fd(i_point,1)] = YD(Pars, n_cosines(i_point,:), n_weight(i_point,1), sigmaT, rho(i_point,1));
    if fd(i_point,1) > 0
        activated(1,i)= i_point;
        i=i+1;
    end
end        

end

function [fd_all] = fdDP1(Pars,rho,drho,n_cosines,n_weight,sigmaT,activated)
%========================================================
% Check convergency rate
%========================================================

[~,n]=size(activated);
fd_act=zeros(n,1);
fd_all=0;

for i=1:n
    i_point=activated(1,i);
    rho(i_point,1)=rho(i_point,1)+drho(i,1);
    [~,fd_act(i,1)] = YD(Pars, n_cosines(i_point,:), n_weight(i_point,1), sigmaT, rho(i_point,1));
    fd_all=fd_all+fd_act(i,1)^2;
end
fd_all=sqrt(fd_all);

end

function [Yd,fd] = YD(Pars, n, w, sigmaT,rho)
%========================================================
%      calculate Yd fd
%========================================================
% global nu0 E0 Kc Ec Ko Eo alpha c0 c1

E0   = Pars(1);
nu0  = Pars(2);
Kc   = Pars(5);
Ec   = Pars(6);
alpha= Pars(7);
Ko   = Pars(8);
Eo   = Pars(9);
c0=(16/3)*(1-nu0^2)/E0;
c1=(32/3)*(1-nu0^2)/(2-nu0)/E0;

M(1:3,1:3,1:3,1:3)=0;
E = eye(3);

if n*sigmaT*n' <= 0
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    M(i,j,k,l) = w*c0*n(i)*n(j)*n(k)*n(l)+w*c1*((0.25*(n(i)*n(k)*E(j,l)+n(i)*n(l)*E(j,k)...
                                   +E(i,k)*n(j)*n(l)+E(i,l)*n(j)*n(k))-n(i)*n(j)*n(k)*n(l)));
                end
            end
        end
    end
    Yd=0.5*Aij_Bij(Aijkl_Bij( M,sigmaT),sigmaT);
    fd=Yd - alpha*trace(sigmaT) - Ko*(1+Eo*rho); 
else
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    M(i,j,k,l) = w*c1*((0.25*(n(i)*n(k)*E(j,l)+n(i)*n(l)*E(j,k)...
                                   +E(i,k)*n(j)*n(l)+E(i,l)*n(j)*n(k))-n(i)*n(j)*n(k)*n(l)));
                end
            end
        end
    end
    Yd=0.5*Aij_Bij(Aijkl_Bij( M,sigmaT),sigmaT);
    fd=Yd - alpha*trace(sigmaT) - Kc*(1+Ec*rho); 
end

end

function [ dg_dsig ] =  Dg_Dsigma(Pars, n, w, sigmaT)
%========================================================
%      calculate Dg_Dsig for single crack set
%========================================================
% global nu0 E0 c0 c1

E0   = Pars(1);
nu0  = Pars(2);

c0=(16/3)*(1-nu0^2)/E0;
c1=(32/3)*(1-nu0^2)/(2-nu0)/E0;


M(1:3,1:3,1:3,1:3)=0;
E = eye(3);


if n*sigmaT*n' <= 0
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    M(i,j,k,l) = w*c0*n(i)*n(j)*n(k)*n(l)+w*c1*((0.25*(n(i)*n(k)*E(j,l)+n(i)*n(l)*E(j,k)...
                                   +E(i,k)*n(j)*n(l)+E(i,l)*n(j)*n(k))-n(i)*n(j)*n(k)*n(l)));
                end
            end
        end
    end
else
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    M(i,j,k,l) = w*c1*((0.25*(n(i)*n(k)*E(j,l)+n(i)*n(l)*E(j,k)...
                                   +E(i,k)*n(j)*n(l)+E(i,l)*n(j)*n(k))-n(i)*n(j)*n(k)*n(l)));
                end
            end
        end
    end
end

dg_dsig=Aijkl_Bij( M,sigmaT);

end

function [g_sig_sum] = Dg_Dsig_Sum(Pars,ddrho,activated,n_cosines,n_weight,sigmaT)
%========================================================
%      calculate Dg_Dsig for all activated crack sets
%========================================================
% global E0 nu0 c0 c1

E0=Pars(1);
nu0=Pars(2);
c0=(16/3)*(1-nu0^2)/E0;
c1=(32/3)*(1-nu0^2)/(2-nu0)/E0;
E = eye(3);

[~,n]=size(activated);

NN(1:3,1:3,1:3,1:3) = 0;
Tri(1:3,1:3,1:3,1:3) = 0;

for ii=1:n
    i_point=activated(1,ii);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    if n_cosines(i_point,:)*sigmaT*n_cosines(i_point,:)'<=0
                        NN(i,j,k,l) = NN(i,j,k,l) + ddrho(ii,1)*n_weight(i_point,1)*c0*n_cosines(i_point,i)*n_cosines(i_point,j)...
                                  *n_cosines(i_point,k)*n_cosines(i_point,l);
                    end
                    Tri(i,j,k,l) = Tri(i,j,k,l) + ddrho(ii,1)*n_weight(i_point,1)*c1*(0.25*(n_cosines(i_point,i)*n_cosines(i_point,k)*E(j,l)+...
                                   n_cosines(i_point,i)*n_cosines(i_point,l)*E(j,k)+E(i,k)*n_cosines(i_point,j)*n_cosines(i_point,l)...
                                   +E(i,l)*n_cosines(i_point,j)*n_cosines(i_point,k))-n_cosines(i_point,i)*n_cosines(i_point,j)...
                                  *n_cosines(i_point,k)*n_cosines(i_point,l));
                end
            end
        end
    end
end
      
matP=NN+Tri;

g_sig_sum=Aijkl_Bij( matP,sigmaT);

end

function [ Df_Dsig,Df_Drho,fd ] =  DF_DSIG_Drho(Pars, n, w, sigmaT,rho,drho)
%========================================================
% Calculate Df_Dsig, Df_Drho, fd for single crack set
%========================================================

% global nu0 E0 Kc Ec Ko Eo alpha c0 c1

E0   = Pars(1);
nu0  = Pars(2);
Kc   = Pars(5);
Ec   = Pars(6);
alpha= Pars(7);
Ko   = Pars(8);
Eo   = Pars(9);

c0=(16/3)*(1-nu0^2)/E0;
c1=(32/3)*(1-nu0^2)/(2-nu0)/E0;

Trsig=sigmaT(1,1)+sigmaT(2,2)+sigmaT(3,3);

M(1:3,1:3,1:3,1:3)=0;
E = eye(3);

if n*sigmaT*n' <= 0
    Df_Drho=-Ko*Eo;
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    M(i,j,k,l) = w*c0*n(i)*n(j)*n(k)*n(l)+w*c1*((0.25*(n(i)*n(k)*E(j,l)+n(i)*n(l)*E(j,k)...
                                   +E(i,k)*n(j)*n(l)+E(i,l)*n(j)*n(k))-n(i)*n(j)*n(k)*n(l)));
                end
            end
        end
    end
else
    Df_Drho=-Kc*Ec;
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    M(i,j,k,l) = w*c1*((0.25*(n(i)*n(k)*E(j,l)+n(i)*n(l)*E(j,k)...
                                   +E(i,k)*n(j)*n(l)+E(i,l)*n(j)*n(k))-n(i)*n(j)*n(k)*n(l)));
                end
            end
        end
    end
end

Temp=Aijkl_Bij( M,sigmaT);

Df_Dsig= Temp - alpha*E;

if n*sigmaT*n' <= 0
    fd = 0.5*Aij_Bij( Temp,sigmaT ) - alpha*Trsig - Ko*(1+Eo*(rho+drho));
else
    fd = 0.5*Aij_Bij( Temp,sigmaT ) - alpha*Trsig - Kc*(1+Ec*(rho+drho));
end

end

function [Inc_instrain] = Inc_In_strain(Pars,ddrho,drho,n_cosines,n_weight,sigmaT, dsig, activated)

% global E0 nu0 c0 c1

E0=Pars(1);
nu0=Pars(2);
c0=(16/3)*(1-nu0^2)/E0;
c1=(32/3)*(1-nu0^2)/(2-nu0)/E0;
E = eye(3);

[~,n]=size(activated);

NN1(1:3,1:3,1:3,1:3) = 0;
Tri1(1:3,1:3,1:3,1:3) = 0;

NN2(1:3,1:3,1:3,1:3) = 0;
Tri2(1:3,1:3,1:3,1:3) = 0;


for ii=1:n
    i_point=activated(1,ii);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    if n_cosines(i_point,:)*sigmaT*n_cosines(i_point,:)'<=0
                        NN1(i,j,k,l) = NN1(i,j,k,l) + drho(ii,1)*n_weight(i_point,1)*c0*n_cosines(i_point,i)*n_cosines(i_point,j)...
                                  *n_cosines(i_point,k)*n_cosines(i_point,l);
                        NN2(i,j,k,l) = NN2(i,j,k,l) + ddrho(ii,1)*n_weight(i_point,1)*c0*n_cosines(i_point,i)*n_cosines(i_point,j)...
                                  *n_cosines(i_point,k)*n_cosines(i_point,l);
                    end
                    Tri1(i,j,k,l) = Tri1(i,j,k,l) + drho(ii,1)*n_weight(i_point,1)*c1*(0.25*(n_cosines(i_point,i)*n_cosines(i_point,k)*E(j,l)+...
                                   n_cosines(i_point,i)*n_cosines(i_point,l)*E(j,k)+E(i,k)*n_cosines(i_point,j)*n_cosines(i_point,l)...
                                   +E(i,l)*n_cosines(i_point,j)*n_cosines(i_point,k))-n_cosines(i_point,i)*n_cosines(i_point,j)...
                                  *n_cosines(i_point,k)*n_cosines(i_point,l));
                    Tri2(i,j,k,l) = Tri2(i,j,k,l) + ddrho(ii,1)*n_weight(i_point,1)*c1*(0.25*(n_cosines(i_point,i)*n_cosines(i_point,k)*E(j,l)+...
                                   n_cosines(i_point,i)*n_cosines(i_point,l)*E(j,k)+E(i,k)*n_cosines(i_point,j)*n_cosines(i_point,l)...
                                   +E(i,l)*n_cosines(i_point,j)*n_cosines(i_point,k))-n_cosines(i_point,i)*n_cosines(i_point,j)...
                                  *n_cosines(i_point,k)*n_cosines(i_point,l));
                end
            end
        end
    end
end
      
matP1=NN1+Tri1;
matP2=NN2+Tri2;

Inc_instrain=Aijkl_Bij( matP1,dsig)+Aijkl_Bij( matP2,sigmaT);

end

function [ddepsin, ddrho, dsig] = fd_lam(Pars, rho, drho, sigmaT,n_cosines,n_weight,activated,R)
%========================================================
%      Iteration solving for ddepsin, ddrho, dsig
%========================================================

global nu0 E0 Kc Ec Ko Eo alpha c0 c1

E0   = Pars(1);
nu0  = Pars(2);
Kc   = Pars(5);
Ec   = Pars(6);
alpha= Pars(7);
Ko   = Pars(8);
Eo   = Pars(9);

c0=(16/3)*(1-nu0^2)/E0;
c1=(32/3)*(1-nu0^2)/(2-nu0)/E0;

[matDz, ~] = matDO2(Pars,rho,drho,activated,n_cosines,n_weight,sigmaT);
Temp=Aijkl_Bkl(matDz,R);

[~,n]=size(activated);
Coef=zeros(n,n);
Fv=zeros(n,1);

for i=1:n
    i_point=activated(1,i);
    [Df_Dsig,Df_Drho,fd ]=DF_DSIG_Drho(Pars, n_cosines(i_point,:), n_weight(i_point,1), sigmaT, rho(i_point,1),drho(i,1));
    Fv(i,1)=fd - Aij_Bij(Df_Dsig, Temp);
    for j=1:n
        j_point=activated(1,j);
        [dg_dsig] =  Dg_Dsigma(Pars, n_cosines(j_point,:), n_weight(j_point,1), sigmaT);
        Coef(i,j) = 2*Aij_Bij( Df_Dsig,Aijkl_Bkl(matDz,dg_dsig));
        if i==j
            Coef(i,j)=Coef(i,j)-Df_Drho;
        end
    end
end

ddrho=Coef\Fv;

[g_sig_sum] = Dg_Dsig_Sum(Pars,ddrho,activated,n_cosines,n_weight,sigmaT);

dsig =  -Aijkl_Bkl(matDz,R+2*g_sig_sum);

[Inc_instrain] = Inc_In_strain(Pars,ddrho,drho,n_cosines,n_weight,sigmaT, dsig, activated);

ddepsin = R+Inc_instrain;

end


