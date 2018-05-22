function Simulation_DM
%Version 140429
%Add: separation factor, some description
%Conrrection: saving name, indentation
%Note: When using SigmaPlot for plotting, function 'surf_colomn' will transform the "surf" variable to vector.

for FOLD_Global = 1    
    clear;clc
    warning off
    global c
end
for FOLD_claim = 1
    Stagecut1 = [];
    yH2 = [];
    recH2 = [];
    Stagecut2 = [];  %re什么是stagecut?
    yCO2 = [];
    recCO2 = [];
    Pp1 = [];
    Pp2 = [];
    Factor_PI_H2 = [];  %re膜的选择性系数 or DMM的分离因数？
    Factor_PEO_CO2 = [];
end

x_Feed = 10:2:50;                          %Axis y     re:axisy和x都是21个元素的行向量（但下面就转置了），用于提供不同的进料流量和其中一种膜丝数量（可得膜面积比）
y_Nf1 = [ 15000:250:19750, 19950 ];        %Axis x     re：这两个向量的元素是每次循环要用的进料的摩尔流率和fiber number

    x_axis = y_Nf1';                       
    y_axis = x_Feed';
    save Process\y_axis.mat y_axis
    save Process\x_axis.mat x_axis         %re下面选择流动形式，就是说上面的代码不分流动形势通用

Flow_Pattern = input('    Flow pattern:     1 --> Counter_Counter    2 --> Co_Co    3 --> Counter_Co    4 --> Co_Counter    ');      %Determine flow pattern

for o = 1:numel(x_Feed)                   %Colomn
    for m = 1:numel(y_Nf1)                %Line
  
        Feed = x_Feed(o);                 %Feed
    
        xH2 = 0.75;                       %Feed H2 mol fraction
    
        Nf_DM1 = y_Nf1(m);                %Fiber number of PI   re膜纤维（膜丝）数量？
        Nf_DM2 = 2e+4-Nf_DM1;             %Fiber number of PEO  re两种fiber number加和等于两万
    
        for FOLD_DM_Feed = 1
            DM_Feed = [                 %Feed stream    re这个变量很重要，别的函数中也要用
                    303;                %T--K
                    1000;               %P--kPa
                    Feed;               %F--kmol/h
                    xH2;                %He
                    1-xH2-0.1;          %CO2
                    0.1;                %CH4
                        ];
            Pp_DM1 = 100;               %PI  permeate pressure
            Pp_DM2 = 100;               %PEO permeate pressure
         
            DM_error = 1e-8;            %Error

            R = DM_Process_input(DM_Feed,Pp_DM1,Pp_DM2,Nf_DM1,Nf_DM2);             %Input function return membrane area ratio (R)

            load Process\DM_Config.mat               
        end

        tic
        for FOLD_DM_Calculation = 1
            if Flow_Pattern == 1
                fprintf('====Counter-Counter====\r')
                DM_Profile = DM_Iteration_Counter_Counter(DM_Feed,DM_error);
                for FOLD_Results_saving = 1                
                    mkdir Process\DM_Profiles\Counter-Counter
                    buff_save = [ 'save Process\DM_Profiles\Counter-Counter\DM_Profile_Counter_Counter_Nf_DM1_',num2str(Nf_DM1),'_F_',num2str(Feed),'_xH2_',num2str(xH2),'.mat DM_Profile' ];
                    eval(buff_save); 

                    Stagecut1(o,m) = DM_Profile(1,1)/Feed;
                    yH2(o,m) = DM_Profile(4,1);
                    recH2(o,m) = DM_Profile(4,1)*DM_Profile(1,1)/Feed/xH2;
                    Stagecut2(o,m) = DM_Profile(2,1)/Feed;
                    yCO2(o,m) = DM_Profile(3+c+2,1);
                    recCO2(o,m) = DM_Profile(3+c+2,1)*DM_Profile(2,1)/Feed/DM_Feed(5);
                    %Nf1, Feed, F1, F2, Fr, y1, y2, xr, Pp1, Pp2
                    Pp1(o,m) = [ DM_Profile(5+3*c,end) ];
                    Pp2(o,m) = [ DM_Profile(6+3*c,end) ];
                    %Separation factor of PI(H2) and PEO(CO2)
                    Factor_PI_H2(o,m) = DM_Profile(3+1,1)/DM_Profile(3+2,1)/(DM_Feed(4)/DM_Feed(5));  %H2
                    Factor_PEO_CO2(o,m) = DM_Profile(3+c+2,1)/DM_Profile(3+c+1,1)/(DM_Feed(5)/DM_Feed(4));  %CO2
                    save Process\Counter-counter.mat Stagecut1  yH2 recH2 Stagecut2 yCO2 recCO2 Pp1 Pp2 Factor_PI_H2 Factor_PEO_CO2
                    Progress = [ o,m]
                end
            end
    
            if Flow_Pattern == 2
                fprintf('====Co-Co====\r')
                DM_Profile = DM_Iteration_Co_Co(DM_Feed,DM_error);
                for FOLD_Results_saving = 1                
                    mkdir Process\DM_Profiles\Co-Co
                    buff_save = [ 'save Process\DM_Profiles\Co-Co\DM_Profile_Co_Co_Nf_DM1_',num2str(Nf_DM1),'_F_',num2str(Feed),'_xH2_',num2str(xH2),'.mat DM_Profile' ];
                    eval(buff_save); 
        
                    Stagecut1(o,m) = DM_Profile(1,end)/Feed;
                    yH2(o,m) = DM_Profile(4,end);
                    recH2(o,m) = DM_Profile(4,end)*DM_Profile(1,end)/Feed/xH2;
                    Stagecut2(o,m) = DM_Profile(2,end)/Feed;
                    yCO2(o,m) = DM_Profile(3+c+2,end);
                    recCO2(o,m) = DM_Profile(3+c+2,end)*DM_Profile(2,end)/Feed/DM_Feed(5);
                    %Nf1, Feed, F1, F2, Fr, y1, y2, xr, Pp1, Pp2
                    Pp1(o,m) = [ DM_Profile(5+3*c,1) ];
                    Pp2(o,m) = [ DM_Profile(6+3*c,1) ];       
                    %Separation factor of PI(H2) and PEO(CO2)
                    Factor_PI(o,m) = DM_Profile(3+1,end)/DM_Profile(3+2,end)/(DM_Feed(4)/DM_Feed(5));  %H2
                    Factor_PEO(o,m) = DM_Profile(3+c+2,end)/DM_Profile(3+c+1,end)/(DM_Feed(5)/DM_Feed(4));  %CO2        
                    save Process\Co-co.mat Stagecut1  yH2 recH2 Stagecut2 yCO2 recCO2 Pp1 Pp2  Factor_PI_H2 Factor_PEO_CO2
                    Progress = [ o,m]
                end
            end
    
            if Flow_Pattern == 3
                fprintf('====Counter-Co====\r')
                DM_Profile = DM_Iteration_Counter_Co(DM_Feed,DM_error);
                for FOLD_Results_saving = 1                 
                    mkdir Process\DM_Profiles\Counter-Co
                    buff_save = [ 'save Process\DM_Profiles\Counter-Co\DM_Profile_Counter_Co_Nf_DM1_',num2str(Nf_DM1),'_F_',num2str(Feed),'_xH2_',num2str(xH2),'.mat DM_Profile' ];
                    eval(buff_save); 
               
                    Stagecut1(o,m) = DM_Profile(1,1)/Feed;
                    yH2(o,m) = DM_Profile(4,1);
                    recH2(o,m) = DM_Profile(4,end)*DM_Profile(1,1)/Feed/xH2;
                    Stagecut2(o,m) = DM_Profile(2,end)/Feed;
                    yCO2(o,m) = DM_Profile(3+c+2,end);
                    recCO2(o,m) = DM_Profile(3+c+2,end)*DM_Profile(2,end)/Feed/DM_Feed(5);
                    %Nf1, Feed, F1, F2, Fr, y1, y2, xr, Pp1, Pp2
                    Pp1(o,m) = [ DM_Profile(5+3*c,end) ];
                    Pp2(o,m) = [ DM_Profile(6+3*c,1) ];
                    %Separation factor of PI(H2) and PEO(CO2)
                    Factor_PI(o,m) = DM_Profile(3+1,1)/DM_Profile(3+2,1)/(DM_Feed(4)/DM_Feed(5));  %H2
                    Factor_PEO(o,m) = DM_Profile(3+c+2,end)/DM_Profile(3+c+1,end)/(DM_Feed(5)/DM_Feed(4));  %CO2        
                    save Process\Counter-co.mat Stagecut1  yH2 recH2 Stagecut2 yCO2 recCO2 Pp1 Pp2 Factor_PI_H2 Factor_PEO_CO2
                    Progress = [ o,m]
                end
            end   
    
            if Flow_Pattern == 4
                fprintf('====Co-Counter====\r')
                DM_Profile = DM_Iteration_Co_Counter(DM_Feed,DM_error);
                for FOLD_Results_saving = 1                  
                    mkdir Process\DM_Profiles\Co-Counter
                    buff_save = [ 'save Process\DM_Profiles\Co-Counter\DM_Profile_Co_Counter_Nf_DM1_',num2str(Nf_DM1),'_F_',num2str(Feed),'_xH2_',num2str(xH2),'.mat DM_Profile' ];
                    eval(buff_save); 
  
                    Stagecut1(o,m) = DM_Profile(1,end)/Feed;
                    yH2(o,m) = DM_Profile(4,end);
                    recH2(o,m) = DM_Profile(4,end)*DM_Profile(1,end)/Feed/xH2;
                    Stagecut2(o,m) = DM_Profile(2,1)/Feed;
                    yCO2(o,m) = DM_Profile(3+c+2,1);
                    recCO2(o,m) = DM_Profile(3+c+2,1)*DM_Profile(2,1)/Feed/DM_Feed(5);
                    %Nf1, Feed, F1, F2, Fr, y1, y2, xr, Pp1, Pp2
                    Pp1(o,m) = [ DM_Profile(5+3*c,1) ];
                    Pp2(o,m) = [ DM_Profile(6+3*c,end) ];        
                    %Separation factor of PI(H2) and PEO(CO2)
                    Factor_PI(o,m) = DM_Profile(3+1,end)/DM_Profile(3+2,end)/(DM_Feed(4)/DM_Feed(5));  %H2
                    Factor_PEO(o,m) = DM_Profile(3+c+2,1)/DM_Profile(3+c+1,1)/(DM_Feed(5)/DM_Feed(4));  %CO2        
                    save Process\Co-counter.mat Stagecut1  yH2 recH2 Stagecut2 yCO2 recCO2 Pp1 Pp2 Factor_PI_H2 Factor_PEO_CO2
                    Progress = [ o,m]
                end
            end
        end
        toc
    end
end





