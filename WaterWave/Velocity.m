clear all;
close all;
% Last modified: 9/24/2008 9:00pm
% This program is written to compute the Fourier Transform of a function.
% Set parameters.
L=25;
N=256;
h=2*L/N;
k=pi./L;
l=N.*pi/(2.*L);
tau=7.2e-5;
g=9.8;
dt=10;
nframe=100;
wf=15;
hf=.8;
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
    surf(x1,x2,real(eta));
    axis([-wf wf -wf wf -hf hf]);
%%%%%%% Here we are going to compute the wave velocity. %%%%%%%
    %%% Fixed observation point at row 128. This set observation point
    %%% travels along x-axis at y=0. %%%
    [max_eta,column]=max(eta(128,:));
    max_column(t)=column;
    %compute distance travel along x-direction where y is fixed.
    dist_x(t)=column.*h-L;  
    real_time(t)=T_step;
end


% Polyfit
poly=polyfit(real_time,dist_x,1);
slope=poly(:,1)
intersect=poly(:,2);
vel_poly=slope.*real_time+intersect;
figure, plot(real_time,dist_x,'-',real_time,vel_poly,'*'), xlabel('real time in second'), ylabel('distance along xaxis in meter');
title('Wavelet function with parameters in [m,s]:L=25 N=256 T=7.2e-5 g=9.8 dt=10 nframe=100');
text(6,-4,' \leftarrow slope=-0.9041', 'FontSize',12)
legend('log plot','linear plot');
% Save data and plots
save(['E:\MATH\MatLab Files\Save Data_9_21_2008\vel_outputdata4' '.mat'])

