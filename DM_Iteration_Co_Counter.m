function f = DM_Iteration_Co_Counter(DM_feed,DM_error)


    %infolder_name = 'Process\DM_initial.mat';
    %if exist(infolder_name,'file')
        %load Process\DM_initial.mat 
    %else
        DM_Profile_i = DM_Initial(DM_feed);
    %end
        Fp1_DM_a = DM_Profile_i(1,end);
        Fp2_DM_a = DM_Profile_i(2,1);
        Fr_DM_a  = DM_Profile_i(3,end);
        
    DM_Profile = DM_chase_Co_Counter(DM_Profile_i,DM_feed);

        Fp1_DM_b = DM_Profile(1,end);
        Fp2_DM_b = DM_Profile(2,1);
        Fr_DM_b  = DM_Profile(3,end);

    Iteration_Count = 1;
    while Error_return(abs((Fp1_DM_a-Fp1_DM_b)/Fp1_DM_b),DM_error) ~=1 ...
           || Error_return(abs((Fp2_DM_a-Fp2_DM_b)/Fp2_DM_b),DM_error)  ~=1 ...
           || Error_return(abs((Fr_DM_a-Fr_DM_b)/Fr_DM_b),DM_error)  ~=1 ...

            Fp1_DM_a = Fp1_DM_b;
            Fp2_DM_a = Fp2_DM_b;
            Fr_DM_a  = Fr_DM_b;
            
        DM_Profile = DM_Pressure_revolve_Co_Counter(DM_Profile);
        DM_Profile = DM_chase_Co_Counter(DM_Profile,DM_feed);
        
            Fp1_DM_b = DM_Profile(1,end);
            Fp2_DM_b = DM_Profile(2,1);
            Fr_DM_b  = DM_Profile(3,end);
            
        fprintf('%d, ',Iteration_Count);    
        Iteration_Count = Iteration_Count + 1;
        if Iteration_Count >= 200
            DM_Profile = zeros(6+5*c,DM_stage);
            break            
        end
    end
    fprintf('\r    Coverged after %d iterations\r',Iteration_Count);
f = DM_Profile;






