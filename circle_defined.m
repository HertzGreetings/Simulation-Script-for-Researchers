% c=circle_defined(x,y,x0,y0,r)
function c=circle_defined(x,y,x0,y0,r)
c=double((x-x0).^2+(y-y0).^2<=r^2);
end