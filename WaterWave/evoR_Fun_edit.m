clear all;
close all;
% Last modified: 8/15/2008 1:--pm
% This program is written to compute the Fourier Transform of a function.
% Set parameters.
L=25;
N=256;
h=2*L/N;
k=pi./L;
l=N.*pi/(2.*L);
tau=72;
g=980;
dt=300;

%%% Code required to make movie
nframe=70;
W=moviein(nframe);
wf=6;
hf=.8;
% [x1,x2,eta] = peaks(L)
% axis tight
% set(gca,'nextplot','replacechildren');

% Perform Evolution using Fourier Transform.
for t=1:nframe
    % %Discretize time step.
    T_step=(t-1)./dt;
    % %Compute g^(z1,z2).
    [x1,x2]=meshgrid([-L:h:L-h]);
    %fgen is the initial condition function
    y=fgen(x1,x2);
    %Here the computation of G(m,n):=Gmn
    Gmn=fgen(x1,x2).*exp(i.*l.*((x1+x2)+2.*L));
    Gmn_shift=fftshift(Gmn);
    % Compute fft2 of yexp gives g^:=g_head
    g_head=fft2(Gmn_shift);
    % %Compute the time term:
    % cos(sqrt(sqrt(z1.^2+z2.^2).*(g+T.*(z1.^2+z.^2))).*t)
    [z1,z2]=meshgrid([1:1:N]);
    zz1=((z1-1)*k)-mean(mean(z1-1)*k);
    zz2=((z2-1)*k)-mean(mean(z2-1)*k);
    %Time=exp(-.0001*(zz2.^2 + zz1.^2)*T_step);
    azz=sqrt(zz1.^2 + zz2.^2);
    Time=cos(sqrt(g*azz + tau*azz.^3)*T_step);
    %Time=exp(i*sqrt(sqrt((zz1).^2+(zz2).^2).*(g+tau.*((zz1).^2+(zz2).^2))).*T_step);
    % etao = Time.*YEXP which is simplified for the computation of ifft2
    %etao=YEXP;
    % %Compute g^^, where we define g_double_^=g_double_head
    g_double_head=Time.*g_head;
    % %Inverse Fourier Transform
    [m,n]=meshgrid([1:1:N]);
    B=exp(-i.*h.*l.*(m+n));
    ETA=B.*ifft2(g_double_head);
    eta=ifftshift(ETA);
    surf(x1,x2,real(eta));
    axis([-wf wf -wf wf -hf hf]);
    %contour(x1,x2,X);
    % This following code for movie generation is written by Dr.Wright.
%     makemovie=0;
%     if makemovie
%         flag=(mod(j,2)==0);
%         if flag
%             k=k+1;
%             legend('U','linear superposition','Location','Northwest');
%             drawnow;
%             M(k)=getframe;
%         end
%     else
%         drawnow;
%     end
%     j=j+1;
% end
drawnow;
W(:,t)=getframe(gcf);
end;
%movie(W)

%map=colormap;
%mpgwrite(W, map,'Water Wave.mpg')