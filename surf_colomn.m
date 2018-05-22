function f = surf_colomn(z)

c=size(z);
z2=[];
for i = 1:c(2)
    z2(c(1)*(i-1)+1:c(1)*i,1) = z(:,i);
end
f= z2;
