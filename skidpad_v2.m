%% Vehicle Skidpad Simulator v2
%% Initially developed by Song Young ik, KAIST Zilzu Racing
% First release, 2024-04-25, v2.0

% Ignored roll steer, (3d geometry effect -> on kpi,caster effect to
% camber, roll steer), roll center change, cog height change effect by anti
% geometry, camber change by longitudinal acceleration
% MUST BE IMPLEMENTED IN FUTURE!

% If you change the suspension geometry, go to camber function and change
% relationship between roll angle and wheel camber

% guys please don't make rc upper than cog
% and they get camber angle as positive
% assume right turn

tyre = MagicFormulaTyre('C:\Users\songy\Desktop\SONG\KAIST\Club\Zilzu\Formula\Suspension\2024tyre.tir');
%% Part 0 for suspension & steering geometry
% body origin is (0, 0, ride_height)
% in case of knuckle center is (0, 0, r)

Pressure = 83000; %tyre pressure in Pa
Density = 1.225; %environment air density kg/m3
Area = 1; %frontal area of car m2
C_d = 1; %drag coefficient of car
wr = 0.23241; %r of tyre m

%unfinished below here, left for future implementing suspension geometry
%front

%{

rhf = 0.09; %front ride height
fkn_u = [0.06637 0 0.08];
fkn_l = [0.0321 0.008 -0.08];
str_f = [0.018 0.062 0];

fl1 = [-0.2475 0.125 0.03];
fl2 = [-0.25375 -0.125 0.055];
fu1 = [-0.285 0.125 0.18];
fu2 = [-0.285 -0.125 0.18];

%rear
rhr = 0.09; %rear ride height
rkn_u = [0.06637 0 0.08];
rkn_l = [0.0321 0.008 -0.08];
str_r = [0.018 0.062 0];

rl1 = [-0.265 0.1 0.06];
rl2 = [-0.265 -0.1 0.06];
ru1 = [-0.265 0.1 0.2];
ru2 = [-0.265 0.1 0.2];

%Part 0 preprocess
fkn_u(3) = fkn_u(3)+wr;
fkn_l(3) = fkn_u(3)+wr;
rkn_u(3) = fkn_u(3)+wr;
rkn_l(3) = fkn_u(3)+wr;
str_f(3) = str_f(3)+wr;
str_r(3) = str_r(3)+wr;
fl1(3) = fl1(3)+rhf;
fl2(3) = fl2(3)+rhf;
fu1(3) = fu1(3)+rhf;
fu2(3) = fu2(3)+rhf;
rl1(3) = rl1(3)+rhr;
rl2(3) = rl2(3)+rhr;
ru1(3) = ru1(3)+rhr;
ru2(3) = ru2(3)+rhr;

%}

%% Part 1 for steer geometry
%must be implemented using coord system in future
lrack = 0.442; %length of rack m
larm = 0.08; %length of arm at wheel sys m
lrod = 0.38574; %length of steering rod m
rev_pinion = 0.08799; %length of pinion per rev, or 2pi*r, m/rev 
d = 0.075; % longitudinal distance from rack to wheel center, always positive dont worry, m
rpinion = rev_pinion/(2*pi);


%% Part 2 for mass transfer
r = 9; %radius of skidpad, m
%A = 1.2; % G force
%beta = 0; % car slip angle, variable
kfw = 50; %front wheel rate N/mm
krw = 100; %rear wheel rate N/mm
Ws = 300; %sprung weight kg
wu = 15; %Unsprung mass of 2 wheel system, kg
rcf = 0.023; %rollcenter front height, m
rcr = 0.042; %rollcenter rear height, m
l = 1.5; %wheelbase, m
tf = 1.2; %front track, m
tr = 1.17; %rear track, m
wfr = 0.5; %sprung cog percentage front wfr rear 1-wfr
cog = 0.35; %cog height, m

wuf = wu;
wur = wu;

%% Part 3 for kpi & caster
kpi = 10; %kpi angle, degree
caster = 3; %caster angle, degree

kpi = deg2rad(kpi);
caster = deg2rad(caster);

%% Part 4 for static toe
toe_f = 0; %front toe angle, degree
toe_r = 0; %rear toe angle, degree

toe_f = deg2rad(toe_f);
toe_r = deg2rad(toe_r);


%% Main running area

%variable car rotational angle 'beta', steering wheel angle 'theta', and G-force 'A'
%iterate to find A that matches. Probabily between 0 and 2?

%beta = deg2rad(-2);
%theta = deg2rad(20);
%[a1, n1] = iterating(beta,theta,tyre,Pressure,Density,wr,Area,C_d,lrack,larm,lrod,rpinion,d,r,kfw,krw,Ws,wuf,wur,rcf,rcr,l,tf,tr,wfr,cog,kpi,caster,toe_f,toe_r)

s_space = [-140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140];
b_space = [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10];

n_space = zeros(length(s_space),length(b_space));
a_space = zeros(length(s_space),length(b_space));

for s = 1:length(s_space)
    for b = 1:length(b_space)
        disp([s,b])
        [a1, n1] = iterating(deg2rad(b_space(b)),deg2rad(s_space(s)),tyre,Pressure,Density,wr,Area,C_d,lrack,larm,lrod,rpinion,d,r,kfw,krw,Ws,wuf,wur,rcf,rcr,l,tf,tr,wfr,cog,kpi,caster,toe_f,toe_r);
        n_space(s,b) = n1;
        a_space(s,b) = a1; 
    end
    plot(a_space(s,1:end),n_space(s,1:end),'Color','red')
    text(a_space(s,(length(b_space)+1)/2),n_space(s,(length(b_space)+1)/2),num2str(s_space(s)),'Color','red')
    hold on
end

for b = 1:length(b_space)
    plot(a_space(1:end,b),n_space(1:end,b),'--','Color','blue')
    text(a_space(1,b),n_space(1,b),num2str(b_space(b)),'Color','blue')
    hold on
end
grid on
ylabel('Oversteer Moment (Nm)') 
xlabel('Lateral Acceleration (G)')

%% Iterating function, get beta, theta, env and return A, N

function [acc,n] = iterating(beta,theta,tyre,Pressure,Density,wr,Area,C_d,lrack,larm,lrod,rpinion,d,r,kfw,krw,Ws,wuf,wur,rcf,rcr,l,tf,tr,wfr,cog,kpi,caster,toe_f,toe_r)

    if beta == 0 && theta == 0
        acc = 0;
        n = 0;
        return
    end
    

    A_temp = 1;
    lateral = 100;
    count = 0;

    fl_angle0 = steering(0, lrack, larm, lrod, rpinion, d, tf);
    fr_angle0 = steering(0, lrack, larm, lrod, rpinion, d, tf);
    fl_angle = steering(-theta, lrack, larm, lrod, rpinion, d, tf);
    fr_angle = steering(theta, lrack, larm, lrod, rpinion, d, tf);
    fl_angle = fl_angle - fl_angle0;
    fr_angle = fr_angle0 - fr_angle;

    while abs(lateral)>50 && count < 50
        count = count + 1;
        A = A_temp;

        %return slip angle of fl, fr, rl, rr
        [wheel_angle_fl,wheel_angle_fr,wheel_angle_rl,wheel_angle_rr,vel_l,vel_r] = slipa(A,r,beta,wfr,tf,tr,l,toe_f,toe_r,fl_angle,fr_angle);
        
        [load_fl,load_fr,load_rl,load_rr,roll] = loading(Ws, wr, wuf, wur, kfw, krw, tf, tr, rcf, rcr, l, cog, wfr, A, beta);
        
        camber_fl = camber(roll, 'f', 'l') + cambersteer(kpi,caster,fl_angle);
        camber_fr = camber(roll, 'f', 'r') + cambersteer(kpi,caster,-fr_angle);
        camber_rl = camber(roll, 'r', 'l');
        camber_rr = camber(roll, 'r', 'r');
        
        %sign convention is opposite in magicformula, minus in inside dir
        %positive for braking force (acc sr is negative)
        Velocity = sqrt(abs(A)*9.8*r);
        %Velocity = 15;

        [t_fl,f_fl,mz_fl] = magicformula(tyre, 0, wheel_angle_fl, load_fl, Pressure, camber_fl, Velocity);
        [t_fr,f_fr,mz_fr] = magicformula(tyre, 0, wheel_angle_fr, load_fr, Pressure, camber_fr, Velocity);
        
        drag_i_f = f_fl*sin(fl_angle+toe_f)+f_fr*sin(fr_angle-toe_f);
        lateral_f = +f_fl*cos(fl_angle+toe_f)+f_fr*cos(fr_angle-toe_f)+t_fl*sin(fl_angle+toe_f)+t_fr*sin(fr_angle-toe_f);
        drag_aero = 0.5*Density*Velocity*Velocity*Area*C_d;
        
        sr0 = -0.2;
        sr1 = 0.2;
        sr_t = 0;
        it_count = 0;
        thrust = 10;
        thrust0 = 0;
        thrust1 = 0;
        thrust_t = 0;

        while abs(thrust)>1 %&it_count < 10
            it_count = it_count + 1;
            if it_count == 1
                sr = sr0;
            elseif it_count == 2
                sr = sr1;
            elseif it_count == 3
                sr = sr_t;
            else
                if thrust_t*thrust0 > 0
                    sr0 = sr_t;
                elseif thrust_t*thrust0 < 0
                    sr1 = sr_t;
                end
                sr_t = (sr0+sr1)/2;
                sr = sr_t;
            end

            %sr_r = ((abs(vel_l)*sr+vel_l)*(vel_r/vel_r)-vel_r)/abs(vel_r);
            sr_r=sr*sign(vel_l)*sign(vel_r);

            [t_rl,f_rl,mz_rl] = magicformula(tyre, sr, wheel_angle_rl, load_rl, Pressure, camber_rl, Velocity);
            [t_rr,f_rr,mz_rr] = magicformula(tyre, sr_r, wheel_angle_rr, load_rr, Pressure, camber_rr, Velocity);
            drag_i_r = +f_rl*sin(toe_f)+f_rr*sin(-toe_f);
            thrust = -t_fl*cos(fl_angle+toe_f)-t_fr*cos(fr_angle-toe_f)+drag_i_f+drag_i_r+drag_aero-t_rl*cos(toe_r)-t_rr*cos(toe_r)-(Ws+wuf+wur)*A*9.8*sin(beta);
            %we removed f=ma so if it is minus, need less thrust and if plus,
            %need more thrust
            if it_count == 1
                thrust0 = thrust;
            elseif it_count == 2
                thrust1 = thrust;
            else
                thrust_t = thrust;
            end

        end
        
        lateral_r = f_rl*cos(toe_r)+f_rr*cos(-toe_r)+t_rl*sin(toe_r)+t_rr*sin(-toe_r);
        lateral = lateral_f + lateral_r - (Ws+wuf+wur)*A*9.8*cos(beta);

        A_temp = (0.8*A + 0.2*(lateral_f+lateral_r)/((Ws+wuf+wur)*9.8*cos(beta)));
        % relaxation factor 0.2

        if count >= 50
            disp([rad2deg(beta),rad2deg(theta), count])
        end

    end

    acc = A;
    n = (wfr*l)*lateral_f - (1-wfr)*l*lateral_r+ (tf/2)*(t_fl*cos(fl_angle+toe_f)-t_fr*cos(fr_angle-toe_f))+ (tr/2)*(t_rl*cos(toe_r)-t_rr*cos(toe_r)) +mz_fl+mz_fr+mz_rl+mz_rr;
    %Positive is -z axis moment, which leads to oversteer

%[FX, FY, MZ, MY, MX] = magicformula(tyre, slipratio, slipangle, FZ, Pressure, Camber, Velocity)
end

%% Slip angle by r

function [afl,afr,arl,arr,vl,vr] = slipa(A,r,beta,wfr,tf,tr,l,toe_f,toe_r,fl_a,fr_a)
    v = sqrt(abs(A)*9.8*r);
    %v = 15;
    omega = sign(A)*v/r;
    %omega = A/v;

    vy = v*sin(beta);
    vx = v*cos(beta);
    %if alpha is positive right, negative left 
    afl = atan((vy+omega*wfr*l)/(vx+omega*tf/2)) - toe_f - fl_a;
    afr = atan((vy+omega*wfr*l)/(vx-omega*tf/2)) + toe_f - fr_a;
    arl = atan((vy-omega*(1-wfr)*l)/(vx+omega*tr/2)) - toe_r;
    arr = atan((vy-omega*(1-wfr)*l)/(vx-omega*tr/2)) + toe_r;

    vl = (vx+omega*tr/2);
    vr = (vx-omega*tr/2);

    %afl = -afl;
    %afr = -afr;
    %arl = -arl;
    %arr = -arr;
    %magicformula use convention of leftward turn as positive, so we need
    %to reverse
end

%% Loading distribution
function [fl,fr,rl,rr, roll] = loading(Ws, wr, wuf, wur, kfw, krw, tf, tr, rcf, rcr, l, cog, wfr, A, beta)

    Ay = A*cos(beta);
    Ax = A*sin(beta);

    Ws=Ws*9.8;
    wuf=wuf*9.8;
    wur=wur*9.8;
    Kf = kfw*1000*tf*tf/2; %front roll stiffness Nm/rad
    Kr = krw*1000*tr*tr/2; %rear roll stiffness
    h2 = abs( (rcr-rcf)*wfr*l-l*cog+l*rcf )/sqrt((rcr-rcf)^2+l^2);
    
    roll = -Ay * (-Ws*h2) / (Kf+Kr-Ws*h2);
    
    fmt = Ay*Ws/tf*( h2*(Kf-(1-wfr)*Ws*h2)/(Kf+Kr-Ws*h2) + (1-wfr)*rcf ) + wuf/tf*wr; %front mass transfer in N
    
    frt = Ay*Ws/tr*( h2*(Kf-wfr*Ws*h2)/(Kf+Kr-Ws*h2) + wfr*rcr ) + wur/tr*wr; %rear mass transfer in N
    
    longt = cog/l * Ax * (Ws+wuf+wur); %longitudinal mass transfer
    
    fl = wuf/2 + Ws*(1-wfr)/2 + fmt/2;
    fl = fl - longt;
    
    fr = wuf/2 + Ws*(1-wfr)/2 - fmt/2;
    fr = fr - longt;
    
    rl = wur/2 + Ws*wfr/2 + frt/2;
    rl = rl + longt;
    
    rr = wur/2 + Ws*wfr/2 - frt/2;
    rr = rr + longt;

    if rr <= 0
        rr = 1;
    end
    if rl <= 0
        rl = 1;
    end
    if fl <= 0
        fl = 1;
    end
    if fr <= 0
        fr = 1;
    end
end

%% Camber from roll angle, temporary function to be replaced

function c = camber(roll, fr, lr)
    %magicformula use positive value as inner inclined, what we usually
    %call negative, so kinda inversed
    if lr == 'l'
        a=1;
    elseif lr == 'r'
        a=-1;
    end

    if fr == 'f'
        c = 0.62*a*rad2deg(roll)-1;
    elseif fr == 'r'
        c = 0.72*a*rad2deg(roll)-0.89;
    end
    c=deg2rad(c);
end

%% steering angle effect to wheel angle
function b = steering(angle, lrack, larm, lrod, rpinion, d, tw)
    p = rpinion*angle;
    l1 = (tw-lrack)/2 - p;
    l2 = sqrt(l1*l1+d*d);
    b = pi/2 - atan(d/l1) - acos((larm*larm+l2*l2-lrod*lrod)/(2*larm*l2));
end

%% cambersteer, input wheel angle, output camber change in radian
function result = cambersteer(kpi, caster, steer)
    result = acos(sin(kpi)*cos(steer)) + kpi + acos(sin(caster)*sin(steer)) - pi;
    %due to magicformula inverse issue
end
