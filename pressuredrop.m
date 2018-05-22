%======压降计算=====
function f=pressuredrop(DM,Pp0,Fp1,v,distance)                               
global c 
load Process\DM_Config.mat
R=8.314;                 
Pp0=Pp0*1000;               %转为国际标准单位
Fp1=Fp1*1000/3600;         %mol/s

if DM == 1
    di = d_DM1(2);
    Pp1=(Pp0^2+256*R*T*Fp1/3.14/Nf_DM1/di^4*distance*v)^0.5/1000;
else
    di = d_DM2(2);
    Pp1=(Pp0^2+256*R*T*Fp1/3.14/Nf_DM2/di^4*distance*v)^0.5/1000;

end


f=Pp1;


