%======ճ�ȼ��㣨n��֣�=====
function f = viscosity(y)
global c 
    load Process\DM_Config.mat
    y  = y;                %���Ũ��
    vpure = v_pure;
    O = zeros(c,c);
%nԪ��ϵճ��ϵ������    
    for i = 1:1:c
        for j = 1:1:c
            Oij = (1+(vpure(i)/vpure(j))^0.5*(MW(j)/MW(i))^0.25)^2/(8*(1+MW(i)/MW(j)))^0.5;
            O(i,j) = Oij;
        end
    end
    vi = zeros(1,c);vj=zeros(c,c);
    for i = 1:c
        vi(i) = y(i)*vpure(i);
        for j = 1:c
            vj(i,j) = y(j)*O(i,j);
        end
    end
    vo=zeros(1,3);
    for i=1:1:c
        vo(i)=vi(i)/sum(vj(i,:));
    end
%�������ճ��
    v = sum(vo);
    f = v;