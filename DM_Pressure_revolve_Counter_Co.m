function f = DM_Pressure_revolve_Counter_Co(Profile_P)
global  c 
    load Process\DM_Config.mat
    %���һ����Co���͵�һ����Counter��
    v1 = viscosity(Profile_P(4:3+c,1));
    v2 = viscosity(Profile_P(4+c:3+2*c,end));
    Profile_P(5+3*c,1) = pressuredrop(1,Pp_DM1,Profile_P(1,1),v1,dz);   
    Profile_P(6+3*c,end) = pressuredrop(2,Pp_DM2,Profile_P(2,end),v2,dz);   

    for j = 2:DM_stage
        
        v1 = viscosity(Profile_P(4:3+c,j));
        v2 = viscosity(Profile_P(4+c:3+2*c,DM_stage+1-j));
        Profile_P(5+3*c,j) = pressuredrop(1,Profile_P(5+3*c,j-1),Profile_P(1,j),v1,dz);  
        Profile_P(6+3*c,DM_stage+1-j) = pressuredrop(2,Profile_P(6+3*c,DM_stage+1-j+1),Profile_P(2,DM_stage+1-j),v2,dz); 


        
    end
        
f = Profile_P;