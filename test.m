function [t,r]=rayleigh_evaporation(R0,Rint,t0,P)

    p0=1e5;%ambient pressure
    rho=1e3;%density of the liquid
    u=1e-6; %viscosity
    S=72e-3; %surface tension
    k=1.4; %polytropic constant
    T0=360; %Initial temperature
    L=2.26e6; %Latent heat of water
    M=0.018; %molecular weight of water

    
    %% calculate vapour pressure
    
    function [p]=p_vapour(T)
        if T<373
            A=8.07131; B=1730.63; C=233.426;
        else
            A=8.14019; B=1810.94; C=244.485;
        end
        p=10^(A-B/(T-273+C));
    end

    %% express Rayleigh equation in a system of 1st order DE
 
    function yy=rayleigh_ode(t,y)
        %y(1): radius; y(2): interface velocity (d.radius/dt); y(3): vapour
        %pressure
        dy(1)=y(2);
        
        %Bubble temperature
        T_b=T0*((R0/y(1))^(3*(k-1)));
        
        %Pressure inside bubble
        p_g0=p0+2*S/R0-p_vapour(T0);
        p_g=p_g0*((R0/y(1))^(3*k));
        p_v=p_vapour(300);
        p_b=p_v+p_g;
        
        %interface acceleration
        dy(2)=((p_b-p0)/rho-1.5*y(2)*y(2)-4*u*y(2)/y(1)-2*S/(rho*y(1)))/y(1);
        
        %pressure change caused by extra evaporation
        dy(3)=P*8.31*T_b/(L*M*4/3*pi*(y(1)^3));
        yy=[dy(1) dy(2) dy(3)]'; 
    end

    [t,y]=ode45(@rayleigh_ode,[0 t0],[Rint,0,p_vapour(T0)]);

    r=y(:,1);
end
