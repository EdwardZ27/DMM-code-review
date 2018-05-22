function Factor
global c
    Flow_Pattern = input('    Flow pattern:     1 --> Counter_Counter    2 --> Co_Co    3 --> Counter_Co    4 --> Co_Counter    ');
    
x_Feed = 10:2:50;  
y_Nf1 = [ 50, 500, 1000:1000:19000, 19500, 19950 ];
Factor_PI = [];
Factor_PEO = [];
for o = 1:numel(x_Feed) 
for m = 1:numel(y_Nf1)

for FOLD_Non = 1  
    Feed = x_Feed(o); 
    
    xH2 = 0.75;
    
    Nf_DM1 = y_Nf1(m);   
    Nf_DM2 = 2e+4-Nf_DM1;
    
    for FOLD_Basic = 1
    DM_Feed = [ 
            303;               %T--K
            1000;               %P--kPa
            Feed;            %F--kmol/h
            xH2;             %He
            1-xH2-0.1;             %CO2
            0.1;     %N2
                   ];
    Pp_DM1 = 100;
    Pp_DM2 = 100;
         
    DM_error = 1e-8;
    end
end

R = DM_Process_input(DM_Feed,Pp_DM1,Pp_DM2,Nf_DM1,Nf_DM2);
load Process\DM_Config.mat 

if Flow_Pattern == 1
    buff_load = [ 'load Process\DM_Profiles\Counter-Counter\DM_Profile_Counter_Counter_Nf_DM1_',num2str(Nf_DM1),'_F_',num2str(Feed),'_xH2_',num2str(xH2),'.mat DM_Profile' ];
    eval(buff_load);
    Factor_PI(o,m) = DM_Profile(3+1,1)/DM_Profile(3+2,1)/(DM_Feed(4)/DM_Feed(5));  %H2
    Factor_PEO(o,m) = DM_Profile(3+c+2,1)/DM_Profile(3+c+1,1)/(DM_Feed(5)/DM_Feed(4));  %CO2
    save Process\Counter-counter_factor.mat Factor_PI Factor_PEO
end
if Flow_Pattern == 2
    buff_load = [ 'save Process\DM_Profiles\Co-Co\DM_Profile_Co_Co_R_',num2str(A_DM1/A_DM2),'_F_',num2str(Feed),'_xH2_',num2str(xH2),'.mat DM_Profile' ];
    eval(buff_load);
    Factor_PI(o,m) = DM_Profile(3+1,end)/DM_Profile(3+2,end)/(DM_Feed(4)/DM_Feed(5));  %H2
    Factor_PEO(o,m) = DM_Profile(3+c+2,end)/DM_Profile(3+c+1,end)/(DM_Feed(5)/DM_Feed(4));  %CO2
    save Process\Co-co_factor.mat Factor_PI Factor_PEO
end
if Flow_Pattern == 3
    buff_load = [ 'load Process\DM_Profiles\Counter-Co\DM_Profile_Counter_Co_R_',num2str(A_DM1/A_DM2),'_F_',num2str(Feed),'_xH2_',num2str(xH2),'.mat DM_Profile' ];
    eval(buff_load);
    Factor_PI(o,m) = DM_Profile(3+1,1)/DM_Profile(3+2,1)/(DM_Feed(4)/DM_Feed(5));  %H2
    Factor_PEO(o,m) = DM_Profile(3+c+2,end)/DM_Profile(3+c+1,end)/(DM_Feed(5)/DM_Feed(4));  %CO2
    save Process\Counter-co_factor.mat Factor_PI Factor_PEO
end
if Flow_Pattern == 4
    buff_load = [ 'load Process\DM_Profiles\Co-Counter\DM_Profile_Co_Counter_R_',num2str(A_DM1/A_DM2),'_F_',num2str(Feed),'_xH2_',num2str(xH2),'.mat DM_Profile' ];
    eval(buff_load);
    Factor_PI(o,m) = DM_Profile(3+1,end)/DM_Profile(3+2,end)/(DM_Feed(4)/DM_Feed(5));  %H2
    Factor_PEO(o,m) = DM_Profile(3+c+2,1)/DM_Profile(3+c+1,1)/(DM_Feed(5)/DM_Feed(4));  %CO2
    save Process\Co-counter_factor.mat Factor_PI Factor_PEO
end

end
end

clear;