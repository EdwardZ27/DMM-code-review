function f = DM_Initial(DM_feed)
%流股向量为: T   K
%           P   kPa
%           F   kmol/h
%           mol frac

    global c 
    load Process\DM_Config.mat
    Dual_Profile_Initial = zeros(6+5*c,DM_stage);
    Dual_Profile_Initial(:,1) = Dual_Initial_Stage_Cal(DM_feed);
    
    for j = 2:DM_stage
        stage_feed = [ DM_feed(1:2);
                       Dual_Profile_Initial(3,j-1);
                       Dual_Profile_Initial(4+2*c:3+3*c,j-1); ];
        Dual_Profile_Initial(:,j) = Dual_Initial_Stage_Cal(stage_feed);
    end
    
f = Dual_Profile_Initial;



    
    
    
    
function f = Dual_Initial_Stage_Cal(stage_feed)
    
%流股向量为：T   K
%           P   kPa
%           F   kmol/h
%           mol frac
    global c
    
    Qi = zeros(c,2);
    load Process\DM_Config.mat

    xrj0 = stage_feed(4:end);
    yi_DM1 = y_cross(xrj0,Ji_DM1,Pf_DM,Pp_DM1*2);
    yi_DM2 = y_cross(xrj0,Ji_DM2,Pf_DM,Pp_DM2*2);
    y1 = sum(yi_DM1);
    y2 = sum(yi_DM2);
    
    yi_DM1 = yi_DM1/y1;
    yi_DM2 = yi_DM2/y2;
    
    Qi_DM1 = zeros(c,1);Qi_DM2 = zeros(c,1);
    for i = 1:c
        Qi_DM1(i) = Ji_DM1(i)/3/829561.4*dA_DM1*(Pf_DM*xrj0(i)-Pp_DM1*yi_DM1(i));
        Qi_DM2(i) = Ji_DM2(i)/3/829561.4*dA_DM2*(Pf_DM*xrj0(i)-Pp_DM2*yi_DM2(i));
    end
    Q_DM1 = sum(Qi_DM1);
    Q_DM2 = sum(Qi_DM2);
    
    Fp_DM1 = Q_DM1;
    Fp_DM2 = Q_DM2;
    
    Fr = stage_feed(3,1)-Fp_DM1-Fp_DM2;
    
    xrj1 = (stage_feed(3,1)*xrj0-Qi_DM1-Qi_DM2)/Fr;
    
f = [
        Fp_DM1;
        Fp_DM2;
        Fr;
        yi_DM1;
        yi_DM2;
        xrj1;
        Pf_DM;
        Pp_DM1;
        Pp_DM2;
        Ji_DM1;
        Ji_DM2;
                ];











    