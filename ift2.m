% 2D fast inverse fourier transform
function g=ift2(G, delta_f) %function g=ift(G,delta_f)
N=size(G,1);
g=ifftshift(ifft2(ifftshift(G)))*(N*delta_f)^2;
end