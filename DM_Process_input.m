function f = DM_Process_input(DM_Feed,Pp1,Pp2,Nf_DM1,Nf_DM2)

global c
for FOLD_Input_data = 1
    c = 3;
    Ji_DM1 = [ 520;  195;   7   ];   %HM,GPU 
    Ji_DM2 = [ 195; 1580; 130   ];   %CO2M,GPU 
    L = 1;
    d_DM1 = [ 500e-6; 250e-6 ];
    d_DM2 = [ 500e-6; 250e-6 ];
    
    %Nf_DM1 = ceil(A_DM1/pi/L/d_DM1(1));
    %Nf_DM2 = ceil(A_DM2/pi/L/d_DM2(1));
    Nf_DM1 = Nf_DM1;
    Nf_DM2 = Nf_DM2;
    
    A_DM1 = pi*L*d_DM1(1)*Nf_DM1;    %re DMM中其中一种膜的面积=pi*膜长度*膜直径*膜丝数量
    A_DM2 = pi*L*d_DM2(1)*Nf_DM2;
    
    Pf_DM = DM_Feed(2);
    Pp_DM1 = Pp1;
    Pp_DM2 = Pp2;

    DM_stage = 100;           %re是分多少单元格？？
    dA_DM1 = A_DM1/DM_stage;
    dA_DM2 = A_DM2/DM_stage;
    
    dz = L/DM_stage;
end

for FOLD_Pressure_cal = 1
    T    = DM_Feed(1);                                                         %实际温度，K
    MW   = [ 2; 44; 16 ];                                          %氢气、CO2、CH4的分子量
    v_pure = [ 8.907e-006;1.489e-005;1.144e-005 ];
end
    fprintf('    PI/PEO (M1/M2) area ratio: %.4f\r',A_DM1/A_DM2);
    fprintf('    ');
    mkdir Process\DM_Configs
    mkdir Process\DM_Profiles
    save Process\DM_Config.mat
    
    buff_save = [ 'save Process\DM_Configs\DM_Config_R_',num2str(A_DM1/A_DM2),'.mat' ];
    eval(buff_save);
    
f = A_DM1/A_DM2;
    