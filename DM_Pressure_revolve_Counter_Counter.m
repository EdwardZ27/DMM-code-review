function f = DM_Pressure_revolve_Counter_Counter(Profile_P)
global  c 
    load Process\DM_Config.mat
    %最后一级（Co）和第一级（Counter）
    v1 = viscosity(Profile_P(4:3+c,1));
    v2 = viscosity(Profile_P(4+c:3+2*c,1));
    Profile_P(5+3*c,1) = pressuredrop(1,Pp_DM1,Profile_P(1,1),v1,dz);   
    Profile_P(6+3*c,1) = pressuredrop(2,Pp_DM2,Profile_P(2,1),v2,dz);   

    for j = 2:DM_stage
        
        v1 = viscosity(Profile_P(4:3+c,j));
        v2 = viscosity(Profile_P(4+c:3+2*c,j));
        Profile_P(5+3*c,j) = pressuredrop(1,Profile_P(5+3*c,j-1),Profile_P(1,j),v1,dz);  
        Profile_P(6+3*c,j) = pressuredrop(2,Profile_P(6+3*c,j-1),Profile_P(2,j),v2,dz); 

    end
        
f = Profile_P;