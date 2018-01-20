clear all;
close all;
% Last modified: 9/1/2008 5:07pm
% This program is written to compute the Fourier Transform of a function.
% Set parameters.
L=30;
N=256;
h=2*L/N;
k=pi./L;
l=N.*pi/(2.*L);
tau=72;
g=980;
dt=100;
nframe=20;
% Perform Evolution using Fourier Transform.
for t=1:nframe
    % Discretize time step. 
    % --Remark-- 
    % Change from "T_step=(t-1)./dt" to "T_step=t./dt" considering log_time
    T_step=t./dt;
    % Compute g(z1,z2).
    [x1,x2]=meshgrid([-L:h:L-h]);
    y=fgen(x1,x2);
    yexp=fgen(x1,x2).*exp(i.*l.*((x1+x2)+2.*L));
    yexp_shift=fftshift(yexp);
    % Compute fft2
    YEXP=fft2(yexp_shift);
    % Compute the time term:
    % cos(sqrt(sqrt(z1.^2+z2.^2).*(g+T.*(z1.^2+z.^2))).*t)
    [z1,z2]=meshgrid([1:1:N]);
    zz1=((z1-1)*k)-mean(mean(z1-1)*k);
    zz2=((z2-1)*k)-mean(mean(z2-1)*k);
    azz=sqrt(zz1.^2 + zz2.^2);
    Time=cos(sqrt(g*azz + tau*azz.^3)*T_step);
    % etao = Time.*YEXP which is simplified for the computation of ifft2
    etao=Time.*YEXP;
    % Inverse Fourier Transform, where there is an extra term, "B", in the
    % calculation.
    [m,n]=meshgrid([1:1:N]);
    B=exp(-i.*h.*l.*(m+n));
    ETA=B.*ifft2(etao);
    eta=ifftshift(ETA);
    %surf(x1,x2,real(eta)), title('evolution of function in real domain');
    %contour(x1,x2,X);
    
    % Compute the maxium eta(x,y,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     Step 1. This estimates the max eta along the row at some time step
%             Max_eta_dx=max(abs(eta));
%     Step 2. This estimate the ultimate max eta along the coloumn of 
%             the maximized row from step 1. 
%             max_eta_dy=max(Max_eta_filt);
%     Step 3. This gives the ultimate max eta at each time step.
%             max_eta_dxdy(t)=max_eta_dy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 
% The above method simplified into the following code using "max(max())".
% For the setup of Y-axis
    max_eta=real(max(max(eta)));
    A(t)=max_eta;
% Log (A)-->
    log_A(t)=real(log(max_eta));
% For the setup of X-axis
    real_time(t)=T_step;
% Log (real_time)-->
    log_time(t)=log(T_step);
end

% Plot before using logorithm
figure(1), plot(real_time,A), xlabel('time'), ylabel('A');
% Logorithm plot & Polyfit
poly=polyfit(log_time,log_A,1)
slope=poly(:,1);
intersect=poly(:,2);
fun_poly=slope.*log_time+intersect;
figure(2), plot(log_time,log_A,'-',log_time,fun_poly,'*'), xlabel('log time'), ylabel('log A');
title('Wavelet function with parameters in units [cm,s]:L=30 N=256 T=72 g=980 dt=100 nframe=250');
text(-2,-1,' \leftarrow slope=-0.8695', 'FontSize',12)
legend('log plot','linear plot');
% Save data and plots
save(['E:\MATH\MatLab Files\Save Data_9_21_2008\outputdata8' '.mat'])




