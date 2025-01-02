% 同态滤波器，D0为高斯滤波的标准差，rL和rH是截止幅度
% 默认参数：D0=80; rL=0.25; rH=2.2; c=2.0; 
% cnliuixn1995@gmail.com
function g = homo_filter(im,D0,s,rL,rH,c)
[size_y,size_x] = size(im);
M = 2*size_x;
N = 2*size_y;
f = log(im+1); % 对图像的灰度取对数
Fp = fft2(f,M,N); % 傅里叶变换之前进行补充0
% 设计滤波器
[v, u] = meshgrid(1:N,1:M);
u = u - floor(M/2); %中心化
v = v - floor(N/2);
D = sqrt(u.^2 + v.^2);
H = 1 - exp(-c*(D./D0^s)); % 高斯高通滤波器
H = (rH - rL) * H + rL;

H= ifftshift(H); %H的反中心化
Gp= Fp .* H;
gp = real(ifft2(Gp)); % 取绝对值容易引起虚部带来的误差
g = gp(1:size_y,1:size_x);
% 取指数
g = exp(g) - 1;
end




