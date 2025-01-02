% Angular spectrum diffraction → Uout = ang_spec_diffr(Uin,lambda,d1,z)；
% Written by Xin Liu: cnliuxin1995@gmail.com
% According to Nyquist sampling theorm: d1 >= sqrt(lambda*z/N)
% Notice: The coordinates ranges of Objective and Observation planes are
% indentical 

function Uout = ang_spec_diffr(Uin,lambda,d1,z)
N = size(Uin,1); % the number of sampling points
% L = N*d1; % diffractive plane range
k = 2*pi/lambda; % propagation wavenumber

df = 1/(N*d1); % space-frequency sampling spacing
[fx,fy] = meshgrid((-N/2:N/2-1)*df); % frequency coordinates
% dk = 2*pi*df; % arc-frequency sampling spacing

H = exp(1i*k*z*sqrt(1-(lambda*fx).^2-(lambda*fy).^2)); % angular spectrum transfer function
Uf = ft2(Uin,d1);
Uout = ift2(Uf .* H , df);
end
