function f = DM_chase_Counter_Co(DM_Profile,DM_feed)

    global c A B

    load Process\DM_Config.mat

for FOLD_Claim_variables = 1
    lo = zeros(c,1);
    Jj_DM1 = zeros(DM_stage,c);
    Jj_DM2 = zeros(DM_stage,c);
    vj_DM1_initial = zeros(DM_stage,c);
    vj_DM2_initial = zeros(DM_stage,c);
    vj_DM1 = zeros(DM_stage,c);
    vj_DM2 = zeros(DM_stage,c);
    Cr_DM1 = zeros(DM_stage,c);       %方便书写的系数
    Cr_DM2 = zeros(DM_stage,c);       %方便书写的系数
    Cp_DM1 = zeros(DM_stage,c);       %方便书写的系数
    Cp_DM2 = zeros(DM_stage,c);       %方便书写的系数
    
    Fpj_DM1 = DM_Profile(1,1:DM_stage)';
    Fpj_DM2 = DM_Profile(2,1:DM_stage)';
    Frj_DM  = DM_Profile(3,1:DM_stage)';         
    Pfj_DM  = DM_Profile(4+3*c,1:DM_stage)';
    Ppj_DM1 = DM_Profile(5+3*c,1:DM_stage)';
    Ppj_DM2 = DM_Profile(6+3*c,1:DM_stage)';
    

    for i = 1:c
        for j = 1:DM_stage
            Jj_DM1(j,i)  = DM_Profile(6+3*c+i,j)/829561.4;
            Jj_DM2(j,i)  = DM_Profile(6+4*c+i,j)/829561.4;
            vj_DM1_initial(j,i)  = DM_Profile(1,j)*DM_Profile(3+i,j);
            vj_DM2_initial(j,i)  = DM_Profile(2,j)*DM_Profile(3+c+i,j);
        end
        lo(i,1) = DM_feed(3)*DM_feed(3+i);
    end
    for i = 1:c
        for j = 1:DM_stage
            Cr_DM1(j,i) = Jj_DM1(j,i)*dA_DM1*Pfj_DM(j)/Frj_DM(j); 
            Cr_DM2(j,i) = Jj_DM2(j,i)*dA_DM2*Pfj_DM(j)/Frj_DM(j); 
            Cp_DM1(j,i) = Jj_DM1(j,i)*dA_DM1*Ppj_DM1(j)/Fpj_DM1(j,1); 
            Cp_DM2(j,i) = Jj_DM2(j,i)*dA_DM2*Ppj_DM2(j)/Fpj_DM2(j,1); 
        end
    end
end

for FOLD_Equation_matrix = 1
%方程形式：Ax=B
%========DM1的计算
    for i = 1:c
        A = zeros(DM_stage,DM_stage);
        B = zeros(DM_stage,1);
        B(1,1) = (-1-1/Cp_DM1(1,i))*lo(i);
        %方程系数第一行
        for j = 1

            A(j,j)   = -(Cp_DM2(j,i)+1)/Cr_DM2(j,i)-1-(Cp_DM2(j,i)+1)/Cr_DM2(j,i)/Cp_DM1(j+1,i)-(Cr_DM1(j+1,i)+1)/Cp_DM1(j+1,i)/Cr_DM2(j+1,i)-1/Cp_DM1(j+1,i)...
                      -(Cr_DM1(j,i)+1)*(Cp_DM2(j,i)+1)/Cp_DM1(j,i)/Cr_DM2(j,i)-1/Cp_DM1(j,i);
            A(j,j+1) = (Cp_DM2(j+1,i)+1)*(Cr_DM1(j+1,i)+1)/Cr_DM2(j+1,i)/Cp_DM1(j+1,i)+1/Cp_DM1(j+1,i);
        end
        for j = 2

            A(j,j-1) = (Cp_DM2(j-1,i)+1)/Cr_DM2(j-1,i)+1/Cr_DM2(j,i)+1+1/Cp_DM1(j+1,i)/Cr_DM2(j,i)+(Cp_DM2(j-1,i)+1)/Cr_DM2(j-1,i)/Cp_DM1(j,i)...
                      +(Cr_DM1(j,i)+1)/Cp_DM1(j,i)/Cr_DM2(j,i)+1/Cp_DM1(j,i);
            A(j,j)   = -(Cp_DM2(j,i)+1)/Cr_DM2(j,i)-1-(Cp_DM2(j,i)+1)/Cr_DM2(j,i)/Cp_DM1(j+1,i)-(Cr_DM1(j+1,i)+1)/Cp_DM1(j+1,i)/Cr_DM2(j+1,i)-1/Cp_DM1(j+1,i)...
                      -(Cr_DM1(j,i)+1)*(Cp_DM2(j,i)+1)/Cp_DM1(j,i)/Cr_DM2(j,i)-1/Cp_DM1(j,i);
            A(j,j+1) = (Cp_DM2(j+1,i)+1)*(Cr_DM1(j+1,i)+1)/Cr_DM2(j+1,i)/Cp_DM1(j+1,i)+1/Cp_DM1(j+1,i);
        end
        for j = 3:DM_stage-1
            A(j,j-2) = -1/Cr_DM2(j-1,i)-1/Cp_DM1(j,i)/Cr_DM2(j-1,i); 
            A(j,j-1) = (Cp_DM2(j-1,i)+1)/Cr_DM2(j-1,i)+1/Cr_DM2(j,i)+1+1/Cp_DM1(j+1,i)/Cr_DM2(j,i)+(Cp_DM2(j-1,i)+1)/Cr_DM2(j-1,i)/Cp_DM1(j,i)...
                      +(Cr_DM1(j,i)+1)/Cp_DM1(j,i)/Cr_DM2(j,i)+1/Cp_DM1(j,i);
            A(j,j)   = -(Cp_DM2(j,i)+1)/Cr_DM2(j,i)-1-(Cp_DM2(j,i)+1)/Cr_DM2(j,i)/Cp_DM1(j+1,i)-(Cr_DM1(j+1,i)+1)/Cp_DM1(j+1,i)/Cr_DM2(j+1,i)-1/Cp_DM1(j+1,i)...
                      -(Cr_DM1(j,i)+1)*(Cp_DM2(j,i)+1)/Cp_DM1(j,i)/Cr_DM2(j,i)-1/Cp_DM1(j,i);
            A(j,j+1) = (Cp_DM2(j+1,i)+1)*(Cr_DM1(j+1,i)+1)/Cr_DM2(j+1,i)/Cp_DM1(j+1,i)+1/Cp_DM1(j+1,i);
        end
        %方程系数最后1行
        for j = DM_stage
            A(j,j-2) = -1/Cr_DM2(j-1,i)-1/Cp_DM1(j,i)/Cr_DM2(j-1,i); 
            A(j,j-1) = (Cp_DM2(j-1,i)+1)/Cr_DM2(j-1,i)+1/Cr_DM2(j,i)+1+(Cp_DM2(j-1,i)+1)/Cr_DM2(j-1,i)/Cp_DM1(j,i)...
                      +(Cr_DM1(j,i)+1)/Cp_DM1(j,i)/Cr_DM2(j,i)+1/Cp_DM1(j,i);
            A(j,j)   = -(Cp_DM2(j,i)+1)/Cr_DM2(j,i)-1 ...
                      -(Cr_DM1(j,i)+1)*(Cp_DM2(j,i)+1)/Cp_DM1(j,i)/Cr_DM2(j,i)-1/Cp_DM1(j,i);

        end
       
        options = optimset( 'display','off' );
        vj_DM2(:,i) = fsolve(@vfun,vj_DM2_initial(:,i),options);
    end

%========DM2的计算
    for i = 1:c
        for j = 1
            vj_DM1(j,i) = -1/Cp_DM1(j,i)*lo(i)...
                         +vj_DM2(j,i)  *((Cp_DM2(j,i)+1)*(Cr_DM1(j,i)+1)/Cr_DM2(j,i)/Cp_DM1(j,i)+1/Cp_DM1(j,i));
        end  
        for j = 2
            vj_DM1(j,i) = vj_DM2(j-1,i)*(-(Cp_DM2(j-1,i)+1)/Cr_DM2(j-1,i)/Cp_DM1(j,i)-(Cr_DM1(j,i)+1)/Cp_DM1(j,i)/Cr_DM2(j,i)-1/Cp_DM1(j,i))...
                         +vj_DM2(j,i)  *((Cp_DM2(j,i)+1)*(Cr_DM1(j,i)+1)/Cr_DM2(j,i)/Cp_DM1(j,i)+1/Cp_DM1(j,i));
        end          
        for j = 3:DM_stage
            vj_DM1(j,i) = vj_DM2(j-2,i)*(1/Cp_DM1(j,i)/Cr_DM2(j-1,i))...
                         +vj_DM2(j-1,i)*(-(Cp_DM2(j-1,i)+1)/Cr_DM2(j-1,i)/Cp_DM1(j,i)-(Cr_DM1(j,i)+1)/Cp_DM1(j,i)/Cr_DM2(j,i)-1/Cp_DM1(j,i))...
                         +vj_DM2(j,i)  *((Cp_DM2(j,i)+1)*(Cr_DM1(j,i)+1)/Cr_DM2(j,i)/Cp_DM1(j,i)+1/Cp_DM1(j,i));
        end        

    end   
end

for FOLD_Revolve_profile = 1
    
    lj_DM  = zeros(DM_stage,c);
    
    xrj_DM = zeros(DM_stage,c);
    yj_DM1 = zeros(DM_stage,c);
    yj_DM2 = zeros(DM_stage,c);
    %===lj====
    for i = 1:c
        %===lj===
        for j = 1
            lj_DM(j,i)  = (Cp_DM2(j,i)+1)/Cr_DM2(j,i)*vj_DM2(j,i);
        end        
        for j = 2:DM_stage
            lj_DM(j,i)  = (Cp_DM2(j,i)+1)/Cr_DM2(j,i)*vj_DM2(j,i)-1/Cr_DM2(j,i)*vj_DM2(j-1,i);
        end

    end

    for j = 1:DM_stage
        Fpj_DM1(j) = sum(vj_DM1(j,:));
        Fpj_DM2(j) = sum(vj_DM2(j,:));
        Frj_DM(j)  = sum(lj_DM(j,:));
    end
    
    for i = 1:c
        for j = 1:DM_stage
            yj_DM1(j,i) = vj_DM1(j,i)/Fpj_DM1(j);
            yj_DM2(j,i) = vj_DM2(j,i)/Fpj_DM2(j);
            
            xrj_DM(j,i)  = lj_DM(j,i)/Frj_DM(j); 
        end
    end

%======revolving profile========    
    %FP1,FP2,Fr
    DM_Profile(1,:)         = Fpj_DM1';
    DM_Profile(2,:)         = Fpj_DM2';
    DM_Profile(3,:)         = Frj_DM';
    %y1,y2,xr
    DM_Profile(4:3+c,:)     = yj_DM1';
    DM_Profile(4+c:3+2*c,:) = yj_DM2';
    DM_Profile(4+2*c:3+3*c,:) = xrj_DM';
    %Pf,Pp1,Pp2
    DM_Profile(4+3*c,:) = Pfj_DM';
    DM_Profile(5+3*c,:) = Ppj_DM1';
    DM_Profile(6+3*c,:) = Ppj_DM2';
    %Ji1,Ji2
    DM_Profile(7+3*c:6+4*c,:) = Jj_DM1'*829561.4;
    DM_Profile(7+4*c:6+5*c,:) = Jj_DM2'*829561.4;
end

f = DM_Profile;

function f = vfun(y)
global A B
f = A*y-B;



