clc;
% clear all;
close all;
disp('Dual SRUKF for BMS--> program start!!')
load("S2_25deg.mat");
% measured_current=FUDS_S2_25deg_80SOC(:,7);
% measured_voltage=FUDS_S2_25deg_80SOC(:,8);
% time=FUDS_S2_25deg_80SOC(:,2);
x_init =[0.91 0 0]';
C1=1005.24554500748;C2=16.5270225350585;R1=0.0497423593776611;
R2=0.0377020407745524;Rs=0.112399209541886;
%C1=1005.24554500748;C2=16.5270225350585;R1=0.0457423593776611;
% R2=0.0277020407745524;Rs=0.112399209541886;
Qact=2;                                  % Cell capacity => 2Ah
openVoltage = 4.2;
tau1=C1*R1;tau2=C2*R2;
theta_init=[Rs R1 tau1 R2 tau2 openVoltage]'; %Parameters
del_soc=[];

Px_init=[0.02 0 0;0 0.5 0;0 0 0.5]; %error covariance matrix
Ptheta_init=[0.01 0 0 0 0 0; 0 0.01 0 0 0 0; 0 0 0.01 0 0 0; ...
             0 0 0 0.01 0 0; 0 0 0 0 0.01 0; 0 0 0 0 0 0.5];

Qx=(10^-8)*eye(3);
Qtheta=(10^-8)*eye(6);
Rx=1;       
Rtheta=10;  
eta = 1;                                                                   % zk = zk_1 - deltaT/Q + SigmaW(1)

%State Function
e=@(tau,deltat)exp(-deltat/tau);
f_x=@(x,theta,Ik,k,deltat)[(x(1)-k*Ik); (e(theta(3),deltat)*x(2) + ...
                            Ik*theta(2)*(1-e(theta(3),deltat))); ...
                            (e(theta(5),deltat)*x(3) + ...
                            Ik*theta(4)*(1-e(theta(5),deltat)))];

%Measurement Function
h_x=@(x,theta,Ik)interp1(soc,ocv,x(1),'linear','extrap') ...
    - x(2) - x(3) - Ik*theta(1);
h_theta=@(x,theta,Ik)theta(6)-x(2)-x(3)-Ik*theta(1);


%Initializing the mean and covariance
x_meas=x_init;
Px_meas=Px_init; Sx_meas=chol(Px_meas);
theta_meas=theta_init;
Ptheta_meas=Ptheta_init; Stheta_meas=chol(Ptheta_meas);

%Variables for plotting
x_meas_LUT2(2)=x_init(1);  %This is for plotting SOC 
open_vol_lut=theta_init(6); 
theta_lut=[];
time_cap=[];
P_cha_max=[];
P_dis_max=[];
delta_soc=0;


lambda_RLS=0.997005;
alpha_x=0.02;
% alpha_theta=0.1;
alpha_theta=0.5;
beta=2;

%----
x_length=length(x_meas);
xsigma_length=2*x_length+1;
theta_length=length(theta_meas);
thetasigma_length=2*theta_length+1;
%----
x_kappa=0;
 theta_kappa=3-x_length;
% theta_kappa=0;
%----
x_lamda=(alpha_x^2)*(x_length+x_kappa)-x_length;
theta_lamda=(alpha_theta^2)*(theta_length+theta_kappa)-theta_length;
%----
x_gamma=sqrt(x_length+x_lamda);
theta_gamma=sqrt(theta_length+theta_lamda);
%----


%BEGIN- WEIGHT COMPUTATION FOR X
% xwm=x_lamda/(x_length+x_lamda);
% xwc=xwm+(1-alpha_x^2+beta);
% for j=1:2*x_length
%     xwm=[xwm 1/(2*(x_length+x_lamda))];
%     xwc=[xwc 1/(2*(x_length+x_lamda))];
% end

xc=x_length+x_lamda;
xwm=[x_lamda/xc 0.5/xc+zeros(1,2*x_length)];           %weights for means
xwc=xwm;
xwc(1)=xwc(1)+(1-alpha_x^2+beta); 

%END- WEIGHT COMPUTATION FOR X

%BEGIN- WEIGHT COMPUTATION FOR THETA
% theta_wm=theta_lamda/(theta_length+theta_lamda);
% theta_wc=theta_wm+(1-alpha_theta^2+beta);
% for j=1:2*theta_length
%     theta_wm=[theta_wm 1/(2*(theta_length+theta_lamda))];
%     theta_wc=[theta_wc 1/(2*(theta_length+theta_lamda))];
% end


thetac=theta_length+theta_lamda;
theta_wm=[theta_lamda/thetac 0.5/thetac+zeros(1,2*theta_length)];           %weights for means
theta_wc=theta_wm;
theta_wc(1)=theta_wc(1)+(1-alpha_theta^2+beta); 
%END- WEIGHT COMPUTATION FOR THETA
check=[];
ax_sigma_hex=[];
tic
for i=3000:13681

Ik=-FUDS_S2_25deg_80SOC(i,7); %PRESENT MEASURED CURRENT(Ik)
Ik_1=-FUDS_S2_25deg_80SOC(i-1,7); %PREVIOUS MEASURED CURRENT(Ik-1)
Vk=FUDS_S2_25deg_80SOC(i,8); %PRESENT VOLTAGE(Vk)
deltat=FUDS_S2_25deg_80SOC(i,2)-FUDS_S2_25deg_80SOC(i-1,2);%STEP TIME


k=eta*(deltat)/(Qact*3600);


%BEGIN--- COMPUTE TIME UPDATE FOR THETA

theta_time=theta_meas;
Stheta_time=(lambda_RLS^-0.5)*Stheta_meas;
%END--- COMPUTE TIME UPDATE FOR THETA

%BEGIN--- COMPUTE TIME UPDATE FOR X

%BEGIN COMPUTE X SIGMA POINTS TIME  
Xsigma_time=x_meas;
for j=1:x_length
    Xsigma_time=[Xsigma_time x_meas+x_gamma*Sx_meas(:,j)];
end
for j=1:x_length
    Xsigma_time=[Xsigma_time x_meas-x_gamma*Sx_meas(:,j)];
end




%END COMPUTE X SIGMA POINTS TIME  

%BEGIN COMPUTE X SIGMA
for j=1:length(Xsigma_time(1,:)) %all rows of sigma
    x_sigma(:,j)=f_x(Xsigma_time(:,j),theta_time,Ik_1,k,deltat);  
end
% ax_sigma_hex=[ax_sigma_hex ',0x'];
% ax_sigma_hex=[ax_sigma_hex num2hex(x_sigma(1,1))];
%END COMPUTE X SIGMA

%BEGIN COMPUTE X TIME
%x_time=(xwm*x_sigma')';
x_time(1,1)=xwm(1)*x_sigma(1,1)+xwm(2)*x_sigma(1,2)+xwm(3)*x_sigma(1,3)+ ...
            xwm(4)*x_sigma(1,4)+xwm(5)*x_sigma(1,5)+xwm(6)*x_sigma(1,6)+ ...
            xwm(7)*x_sigma(1,7);
x_time(2,1)=xwm(1)*x_sigma(2,1)+xwm(2)*x_sigma(2,2)+xwm(3)*x_sigma(2,3)+ ...
            xwm(4)*x_sigma(2,4)+xwm(5)*x_sigma(2,5)+xwm(6)*x_sigma(2,6)+ ...
            xwm(7)*x_sigma(2,7);
x_time(3,1)=xwm(1)*x_sigma(3,1)+xwm(2)*x_sigma(3,2)+xwm(3)*x_sigma(3,3)+ ...
            xwm(4)*x_sigma(3,4)+xwm(5)*x_sigma(3,5)+xwm(6)*x_sigma(3,6)+ ...
            xwm(7)*x_sigma(3,7);
%END COMPUTE X TIME
%BEGIN COMPUTE WEIGHTED DIFF SIGMA
d=(x_sigma-repmat(x_time,1,xsigma_length))';
d1=sqrt(xwc(2))*(d(2:xsigma_length,:));
%END COMPUTE WEIGHTED DIFF SIGMA

%BEGIN QR DECOMPOSITION AND CHOLESKY DOWNDATE
[Q,Sx_time]=qr(d1);
%   d2=[0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;];
%   d1=[d1 d2];
% Sx_time=-traingularization_6d(d1);
Sx_time=Sx_time(1:3,1:3)+chol(Qx);
%Sx_time=Sx_time(1:3,:)+chol(Qx);
%[Sx_time,xerr_msg_time]=cholupdate(Sx_time,(sqrt(abs(xwc(1)))*d(1,:))','-');
% [Sx_time,xerr_msg_time]=mycholupdowndate(Sx_time',(sqrt(abs(xwc(1)))*d(1,:))','-');
[Sx_time,xerr_msg_time]=mycholupdowndate_cordic(Sx_time',(sqrt(abs(xwc(1)))*d(1,:))','-');


if(xerr_msg_time~=0)
    xerr_msg_time;
    break;
end
%END QR DECOMPOSITION AND CHOLESKY DOWNDATE
     
%END--- COMPUTE TIME UPDATE FOR X

%BEGIN- MEASUREMENT UPDATE FOR THETA

%BEGIN THETA SIGMA MEAS
thetaSigma_meas=theta_time;
for j=1:theta_length
    thetaSigma_meas=[thetaSigma_meas theta_time+theta_gamma*Stheta_time(:,j)];
end
for j=1:theta_length
    thetaSigma_meas=[thetaSigma_meas theta_time-theta_gamma*Stheta_time(:,j)];
end
%END THETA SIGMA MEAS

%BEGIN D SIGMA MEAS
for j=1:thetasigma_length %all rows of sigma
    D_sigma(:,j)=h_theta(x_time,thetaSigma_meas(:,j),Ik);  
end
%END D SIGMA END
% D_mean=(theta_wm*D_sigma')';
D_mean=theta_wm(1)*D_sigma(:,1)+theta_wm(2)*D_sigma(:,2)+ ...
       theta_wm(3)*D_sigma(:,3)+theta_wm(4)*D_sigma(:,4)+ ...
       theta_wm(5)*D_sigma(:,5)+theta_wm(6)*D_sigma(:,6)+ ...
       theta_wm(7)*D_sigma(:,7)+theta_wm(8)*D_sigma(:,8)+ ...
       theta_wm(9)*D_sigma(:,9)+theta_wm(10)*D_sigma(:,10)+ ...
       theta_wm(11)*D_sigma(:,11)+theta_wm(12)*D_sigma(:,12)+ ...
       theta_wm(13)*D_sigma(:,13);
d=(D_sigma-repmat(D_mean,1,thetasigma_length))';
d1=sqrt(theta_wc(2))*d(2:thetasigma_length,:);
[Q,Sd]=qr(sqrt(theta_wc(2))*d(2:thetasigma_length,:));
% Sd = vect_10d(d1(1),d1(2),d1(3),d1(4),d1(5),d1(6),d1(7),d1(8),d1(9),d1(10));
Sd=Sd(1,:)+sqrt(Rtheta);
% [Sd,derr_msg]=cholupdate(Sd,(sqrt(abs(dwc(1)))*d(1,:))','-');
% [Sd,derr_msg]=mycholupdowndate(Sd',(sqrt(abs(dwc(1)))*d(1,:))','-');
%check=[check Sd];


[Sd,derr_msg]=mycholupdowndate_cordic(Sd',(sqrt(abs(theta_wc(1)))*d(1,:))','-');
if(derr_msg~=0)
    derr_msg;
    break;
end

%BEGIN: COMPUTE CROSS CORRELATION FOR THETA
d1=thetaSigma_meas-repmat(theta_time,1,thetasigma_length);
d2=D_sigma-repmat(D_mean,1,thetasigma_length);
PthetaD=zeros(6,1);
for j=1:thetasigma_length
    PthetaD=PthetaD+theta_wc(j)*d1(:,j)*d2(:,j)';
end
%END: COMPUTE CROSS CORRELATION FOR THETA
if(Sd~=0)
kk_theta=(PthetaD/Sd')/Sd;
else
    kk_theta=zeros(5,1);
end

%BEGIN: THETA MEAN AND COVARIANCE MEASUREMENT UPDATE
theta_meas=theta_time+kk_theta*(Vk-D_mean);
Ud=kk_theta*Sd;
% [Stheta_meas,thetaerr_msg_update]=cholupdate(Stheta_time,Ud,'-');
% [Stheta_meas,thetaerr_msg_update]=mycholupdowndate(Stheta_time',Ud,'-');
[Stheta_meas,thetaerr_msg_update]=mycholupdowndate_cordic(Stheta_time',Ud,'-');
% Sx_meas=Sx_time*Sx_time'-Uy*Uy';
if(thetaerr_msg_update~=0)
    thetaerr_msg_update;
    break;
end
%BEGIN: THETA MEAN AND COVARIANCE MEASUREMENT UPDATE

%END -MEASUREMENT UPDATE FOR THETA

%BEGIN MEASUREMENT UPDATE FOR X
%BEGIN X SIGMA MEAS
Xsigma_meas=x_time;
for j=1:x_length
    Xsigma_meas=[Xsigma_meas x_time+x_gamma*Sx_time(:,j)];
end
for j=1:x_length
    Xsigma_meas=[Xsigma_meas x_time-x_gamma*Sx_time(:,j)];
end
%END X SIGMA MEAS

for j=1:xsigma_length 
    y_sigma(:,j)=h_x(Xsigma_meas(:,j),theta_time,Ik);  
end
% y_mean=(xwm*y_sigma')';
y_mean=xwm(1)*y_sigma(:,1)+xwm(2)*y_sigma(:,2)+xwm(3)*y_sigma(:,3)+xwm(4)*y_sigma(:,4)+xwm(5)*y_sigma(:,5)+xwm(6)*y_sigma(:,6)+xwm(7)*y_sigma(:,7);
d=(y_sigma-repmat(y_mean,1,xsigma_length))';
d1=sqrt(xwc(2))*d(2:xsigma_length,:);
[Q,Sy]=qr(sqrt(xwc(2))*d(2:xsigma_length,:));
% Sy=vect_6d(d1(1),d1(2),d1(3),d1(4),d1(5),d1(6));
Sy=Sy(1,:)+sqrt(Rx);
% [Sy,yerr_msg]=cholupdate(Sy,(sqrt(abs(xwc(1)))*d(1,:))','-');
% [Sy1,yerr_msg]=mycholupdowndate(Sy',(sqrt(abs(xwc(1)))*d(1,:))','-');
[Sy,yerr_msg]=mycholupdowndate_cordic(Sy',(sqrt(abs(xwc(1)))*d(1,:))','-');
% Sy=Sy1(1,1);
if(yerr_msg~=0)
    yerr_msg;
    break;
end

d3=Xsigma_meas-repmat(x_time,1,xsigma_length);
d4=y_sigma-repmat(y_mean,1,xsigma_length);

Pxy=zeros(3,1);
wtt=[];
for j=1:xsigma_length
    %wtt=[wtt xwc(j)*d3(:,j)*d4(:,j)']
    Pxy=Pxy+xwc(j)*d3(:,j)*d4(:,j)';
end

if(Sy~=0)
kk_x=(Pxy/Sy')/Sy;
else
    kk_x=zeros(3,1);
    break;
end

% if i==100
%     gain=gain+1;
% end
% kk=(Pxy/Sy')/Sy;
x_meas=x_time+kk_x*(Vk-y_mean);
Uy=kk_x*Sy;
% [Sx_measa,bxerr_msg_update]=cholupdate(Sx_time,Uy,'-');
% Sx_measa=[Sx_measa; Sx_measa]
%  [Sx_meas,xerr_msg_update]=mycholupdowndate(Sx_time',Uy,'-');
[Sx_meas,xerr_msg_update]=mycholupdowndate_cordic(Sx_time',Uy,'-');
if(xerr_msg_update~=0)
    xerr_msg_update;
    break;
end
%END MEASUREMENT UPDATE FOR X


x_meas_LUT2(i-2998+1)=x_meas(1)+0.06;
%     x_meas_LUT2(i-2998+1)=x_meas(1);
if(i<12000)
    delta_soc=delta_soc-Ik*deltat;
    del_soc=[del_soc, delta_soc]; 
    open_vol_lut=[open_vol_lut; theta_meas(6)];
    time_cap=[time_cap;FUDS_S2_25deg_80SOC(i,2)];
    
end
theta_lut=[theta_lut theta_meas];
I_dis_max=(theta_meas(6)-2.7)/(theta_meas(1)+theta_meas(2)+theta_meas(4));
I_cha_max=(4.2-theta_meas(6))/(theta_meas(1)+theta_meas(2)+theta_meas(4));
P_dis_max=[P_dis_max I_dis_max*2.7];
P_cha_max=[P_cha_max I_cha_max*4.2];
end



figure 
subplot(2,1,1)
t_soc = true_soc;
ocv_soh=interp1(soc,ocv,t_soc,'linear','extrap');
plot(open_vol_lut(2:9001),'LineWidth',1.75)
hold on
plot(ocv_soh(2:9001),'LineWidth',1.75)
title('OCV Estimation Using Dual SRUKF','color','k')
xlabel('Number of Sample','color','b')
ylabel('OCV(v)','color','b')
legend({'Co-Estimated OCV','Reference OCV'})
% xtickangle(45)
axesH = gca;
axesH.XAxis.TickLabelInterpreter = 'latex';
axesH.XAxis.TickLabelFormat      = '\\textbf{%g}';
axesH.YAxis.TickLabelInterpreter = 'latex';
axesH.YAxis.TickLabelFormat      = '\\textbf{%g}';
% axesH.FontSize = 15;
set(gca,'FontName','Verdana','Fontsize',14,'XColor','k','YColor','k')




soc1=0.7565;% at 3000th data vol 3.91???
soc2=interp1(ocv,soc,open_vol_lut(9001),'linear','extrap');
soc2_estim=x_meas_LUT2(9001);
% soc2=x_meas_LUT2(9001)
del_charge=delta_soc/3600;
del_charge_n=del_soc(5001)/3600;
cap=del_charge/(soc2-soc1);
cap_estim=del_charge/(soc2_estim-soc1);
soc3=interp1(ocv,soc,open_vol_lut(5001),'linear','extrap');
cap_nw=del_charge_n/(soc2-soc1);




toc
%  x_meas_LUT3=x_meas_LUT3+0.06;
%  x_meas_LUT2=x_meas_LUT2+0.06;
subplot(2,1,2)
plot((t_soc)*100,'k--','LineWidth',1.75);
% plot(PRC_SOC(2:5634,:)*100,);
hold on
% % plot(x_meas_LUT3*100,'LineWidth',1.75);
hold on
plot(x_meas_LUT2*100,'LineWidth',1.75);
% plot(x_meas_LUT(:,2:9000)*100,'m','LineWidth',1.75);
% plot(x_meas_LUT1(:,2:9000)*100,'r','LineWidth',1.75);

%plot(x_meas_LUT(:,1297:3361)+(0.2939-0.1099),'LineWidth',1.75);
 title('SOC ESTIMATION BASED ON Dual SR-UKF','color','k')
xlabel('Number of Sample','color','b')
ylabel('SOC(%) ','color','b')
legend({'Referene','D-SRUKF ESTIMATED','D-SRUKF Co-ESTIMATED',})
% legend({'Referene','D-SRUKF ESTIMATED',})
% xtickangle(45)
axesH = gca;
% axesH.XAxis.XLim=[0 9];
axesH.XAxis.TickLabelInterpreter = 'latex';
axesH.XAxis.TickLabelFormat      = '\\textbf{%g}';
axesH.YAxis.TickLabelInterpreter = 'latex';
axesH.YAxis.TickLabelFormat      = '\\textbf{%g}';
% axesH.FontSize = 15;
set(gca,'FontName','Verdana','Fontsize',14,'XColor','k','YColor','k')
set(gca, 'XDir')
legend({'Referene','D-SRUKF ESTIMATED','D-SRUKF Co-ESTIMATED',})
Req_power=[123 118 113 117 117 118 124 119 113];
soc_power=[10 20 30 40 50 60 70 80 90];
ocv_power=interp1(soc,ocv,soc_power/100,'linear','extrap');
Iref_dis_max=(ocv_power-2.7)./(Req_power*0.001);
Iref_cha_max=(4.2-ocv_power)./(Req_power*0.001);
Pref_dis_max=Iref_dis_max*2.7;
Pref_cha_max=Iref_cha_max*4.2;
figure
subplot(2,1,1)
% plot(soc_power,Pref_dis_max,':o',t_soc(3:10684)*100,P_dis_max,':.','LineWidth',1.75);
% hold on
plot(soc_power,Pref_cha_max,':o',t_soc(3:10684)*100,P_cha_max,':.','LineWidth',1.75);
title('Estimation of SOP','color','k')
xlabel('SOC(%)','color','b')
ylabel('Power(watt)','color','b')
legend({'Reference','Co-Estimated'},'Location','northwest')
% legend({'Reference discharge power','Estimated discharge power'})
% xtickangle(45)
axesH = gca;
axesH.XAxis.TickLabelInterpreter = 'latex';
axesH.XAxis.TickLabelFormat      = '\\textbf{%g}';
axesH.YAxis.TickLabelInterpreter = 'latex';
axesH.YAxis.TickLabelFormat      = '\\textbf{%g}';
% axesH.FontSize = 15;
set(gca,'FontName','Verdana','Fontsize',14,'XColor','k','YColor','k')
set(gca, 'XDir','reverse')


xyz=[P_cha_max(9320) P_cha_max(7951), P_cha_max(6581) P_cha_max(5215) P_cha_max(3830) P_cha_max(2389) P_cha_max(1003) 0 0];
figure
subplot(2,1,1)
plot(soc_power,Req_power,':o',t_soc(1:10000)*100,(theta_lut(1,1:10000)+theta_lut(2,1:10000)+theta_lut(4,1:10000))*1000,':.','LineWidth',1.75);
title('Estimation of Equivalent Resistance','color','k')
xlabel('SOC(%)','color','b')
ylabel('Req(milli ohm)','color','b')
legend({'Reference','Co-Estimated'})
% xtickangle(45)
axesH = gca;
axesH.XAxis.TickLabelInterpreter = 'latex';
axesH.XAxis.TickLabelFormat      = '\\textbf{%g}';
axesH.YAxis.TickLabelInterpreter = 'latex';
axesH.YAxis.TickLabelFormat      = '\\textbf{%g}';
% axesH.FontSize = 15;
set(gca,'FontName','Verdana','Fontsize',14,'XColor','k','YColor','k')
set(gca, 'XDir','reverse')

cap;
cap_nw;
cap_estim;
cosrukf_mae = mae(x_meas_LUT2(:,2:9001)-(t_soc(2:9001,:) )'); %0.0105
improved_err=((0.1047-0.0105)/0.1047)*100;