function f = y_cross(xr,Jij,Pf_DM,Pp)
global c xr_cross J_cross Pp_cross c Pf
xr_cross = xr;
J_cross = Jij/829561.4;
Pp_cross = Pp;
Pf = Pf_DM;
options = optimset( 'display','off' );
y_initial = xr;
y_results = fsolve(@y_equations,y_initial,options);
y_crossflow = y_results;

f = y_crossflow;

%======y'计算方程组======
function f = y_equations (y)
    global xr_cross J_cross Pp_cross c Pf
    J_c = J_cross;
    O = zeros(c,1);
    S = zeros(c,1);
    for i = 1:c
        O(i,1) = J_c(i,1)*(Pf*xr_cross(i,1)-Pp_cross*y(i,1));
    end
    for j = 1:c
        S(j) = O(j)/sum(O(:,:));
    end
f = y-S;
