% This script was written by XinLiu: cnliuxin1995@gmail.com
% Website:
% https://scholar.google.com.hk/citations?user=V2zRkHAAAAAJ&hl=zh-CN;
% https://www.researchgate.net/profile/Xin-Liu-277;

clc;
close all
clear
H=800;
D=200;
d1=10e-3;
lambda=633e-6; k=2*pi/lambda;
xlist=(-H/2:H/2-1)*d1; ylist=(-H/2:H/2-1)*d1;
[x,y]=meshgrid(xlist,ylist);

df=1/(H*d1);
fxlist=(-H/2:H/2-1)*df; fylist=(-H/2:H/2-1)*df; 
[fx1,fy1]=meshgrid(fxlist,fylist);

Cn = randn(H) + 1i*randn(H);
RG = 12;
Random_Phase = exp(1i*angle(ft2(Cn.*exp(-(fx1.^2+fy1.^2)/RG^2),1)));

%%
Obj1=zeros(H);
imageA = imread('5.bmp');
imageA = double(imageA); % 灰度图像
imageA=imresize(imageA,[D D],'nearest');
[picpx,picpy]=size(imageA);
picpxnum=round((H-picpx)/2);picpynum=round((H-picpy)/2);
Obj1(picpxnum:picpxnum+picpx-1,picpynum:picpynum+picpy-1)=rot90(imageA,-2);

coh = 15;
pA = exp(-(fx1.^2+fy1.^2)/coh^2);
%%
I3 = zeros(H);
for ii=1:200
Rn = sqrt(pA).*(1/sqrt(2)*(randn(H)+1i*randn(H))); % image A
E1 = ft2(Rn,df);

u1 = E1 .* Obj1;
u = 500;
u2 = ang_spec_diffr(u1,lambda,d1,u);
v=200;
u3 = ang_spec_diffr(u2.*Random_Phase,lambda,d1,v);
I3 = I3 + abs(u3).^2;
ii
end
I3 = I3-min(min(I3));
I3 = I3./max(max(I3));
I3_filter1 = homo_filter(I3,120,4,4,2,0); % 调用同态滤波函数
% -------------------------------------------------------------------------
I3 = zeros(H);

%%
I3_filter = I3_filter1;
C1 = ft2(I3_filter,1);
C2 = abs(ift2(C1.*conj(C1),1));
C2 = C2 - min(min(C2));
C2 = C2./max(max(C2));
load cameracolor.mat
W = nthroot(abs(fftshift(fft2(fftshift(C2.*exp(-(x.^2+y.^2)/5^2))))),2).*circle_defined(x,y,0,0,0.8);
figure(2)
pcolor(xlist,ylist,W);
shading interp
axis square;axis off; colormap('parula')
axis([-0.6 0.6 -0.6 0.6])
clim([0 80])
colorbar

figure(3)
pcolor(xlist,ylist,I3_filter);
shading interp
axis square;axis off; colormap('parula')
colorbar

%% G-S althorithm
coh_four_meas = W;
gs_guss = abs(randn(H) + 1i* randn(H)); % 不用随机数更好
% gs_guss = exp(-x.^2-y.^2); % 不用随机数更好
gs_initial = gs_guss;
beta_list = fliplr((0:0.05:2));
for mm= 1:length(beta_list)
    beta = beta_list(mm);
    for ii=1:20
        coh_four = ft2(gs_initial,df);
        coh_four_pha = angle(coh_four);
        coh_four = coh_four_meas.*exp(1i*coh_four_pha); % 将测量的幅值替换计算的幅值，相位保持不变。
        gs_calcu = (ift2(coh_four,d1));
        gs_initial = gs_calcu.*(imag(gs_calcu)==0&real(gs_calcu)>=0) + (gs_initial - beta*gs_calcu).*(imag(gs_calcu)~=0|real(gs_calcu)<0);
    end
    mm
end
%
gs_initial2 = gs_calcu;
for jj=1:200
    coh_four2 = ft2(gs_initial2,df);
    coh_four_pha2 = angle(coh_four2);
    coh_four2 = coh_four_meas.*exp(1i*coh_four_pha2); % 将测量的幅值替换计算的幅值，相位保持不变。
    gs_calcu2 = (ift2(coh_four2,d1));
    gs_initial2 = gs_calcu2.*(imag(gs_calcu2)==0&real(gs_calcu2)>=0) + 0.*(imag(gs_calcu2)~=0|real(gs_calcu2)<0);
    jj
end


figure(4)
imagesc(fxlist,fylist,fliplr(fftshift(abs(gs_calcu2),1)));
axis xy; axis square; axis off
colormap('hot')

figure(5)
imagesc(fxlist,fylist,fliplr((abs(gs_calcu2))));
axis xy; axis square; axis off




% This script was written by XinLiu: cnliuxin1995@gmail.com  