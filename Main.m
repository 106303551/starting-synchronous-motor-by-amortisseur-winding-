clc
clear

%電機參數設定
r_stat_in=15;%stator內圈半徑
r_stat_out=17;%stator外圈半徑
width_rot=2;%rotor細部寬度
length_rot=14;%rotor細部長度
r=1.5;%rotor半圓半徑
x=r*cosd(0:1:180);%rotor半圓的x軸分量
y=r*sind(0:1:180);%rotor半圓的y軸分量
r_rot=length_rot/2+r;%rotor半徑
theta=0;%rotor初始角度

%stator設定
Bs=12;%stator磁場純量
ang_freq=2*pi;%角頻率(設定定部磁場轉一圈一秒)
dt=0.01%時間片段長度
t=0:dt:25;%總時間長度

%設定stator座標
xy_stat_in=[r_stat_in*cosd(0:4:360);r_stat_in*sind(0:4:360)];%stator內圈xy座標
xy_stat_out=[r_stat_out*cosd(0:4:360);r_stat_out*sind(0:4:360)];%stator外圈xy座標

%設定rotor的xy座標
xy_rot=[[width_rot/2,width_rot/2,width_rot/2+0.5,x,-width_rot/2,-width_rot/2,-width_rot/2,-x,width_rot/2,width_rot/2];
    [0,length_rot/2,length_rot/2,length_rot/2+y,length_rot/2,0,-length_rot/2,-length_rot/2-y,-length_rot/2,0]];
[pol_rot(1,:),pol_rot(2,:)]=cart2pol(xy_rot(1,:),xy_rot(2,:));%將rotor的xy軸座標轉成極座標(用於方便計算)
Bs_x=Bs*cos(ang_freq*t+pi/2);%stator磁場的x分量
Bs_y=Bs*sin(ang_freq*t+pi/2);%stator磁場的y分量
Bs_xy=[Bs_x;Bs_y];%stator磁場xy分量矩陣


%參數設定
R_winding=2*10^(-3);%阻尼線繞組電阻
G=1;%幾何係數
mu=1.2566*10^(-6);%導磁率(銅)
k=r_rot*r_rot*pi*G/mu;%電機結構因數 k=A(r^2*pi)*G/mu
moment_inertia=10^10;%轉動慣量

%設定用於儲存隨時間改變的物理量的矩陣

e_ind=zeros(1,size(t,2));%感應電壓
i_ind=zeros(1,size(t,2));%感應電流
Bw=zeros(1,size(t,2));%轉部磁場
torque_ind=zeros(1,size(t,2));%感應轉矩
ang_accel=zeros(1,size(t,2));%角速度
ang_vel_rot=zeros(1,size(t,2));%rotor轉速

%計算會隨時間改變的物理量(以dt=0.01為計算的時間間隔)
for i=0:1:size(t,2)-1 %dt的次數
    %i+1為目前當下的物理量 i+2為下個dt的物理量
    [vel_x,vel_y]=pol2cart(theta,r_rot*(ang_freq-ang_vel_rot(i+1)));%計算條棒相對磁場在xandy上的速度分量
    e_ind(i+1)=(vel_x*Bs_xy(2,i+1)-vel_y*Bs_xy(1,i+1))*length_rot;%感應電壓=(V外積B)*l(1*2矩陣和1*2矩陣外積即為交叉相乘)
    i_ind(i+1)=e_ind(i+1)/R_winding;%計算感應電流
    Bw(i+1)=mu*i_ind(i+1)/G;%計算繞組磁場的純量
    torque_ind(i+1)=k*((Bw(i+1)*cos(theta))*Bs_xy(2,i+1)-(Bw(i+1)*sin(theta))*Bs_xy(1,i+1));%感應轉矩=k*(繞組磁場外積stat磁場)(交叉相乘)
    ang_accel(i+1)=torque_ind(i+1)/moment_inertia;%角加速度
    
    if i<size(t,2)-1 %為避免多計算一個dt的轉速 用此做限制(確保t and ang_vel_rot矩陣長度一樣)
        ang_vel_rot(i+2)=ang_vel_rot(i+1)+ang_accel(i+1)*dt;%計算隨著dt改變的轉速(下個dt的轉速為原始轉速+(角加速度*dt))
    end
    
    theta_plus=ang_vel_rot(i+1)*dt;%計算這個dt改變的角度(視在這個dt期間轉速為常數)
    theta=theta+theta_plus;%計算下個dt時的角度(此dt時的角度+此dt改變的角度)
    [xy_rot(1,:),xy_rot(2,:)]=pol2cart(pol_rot(1,:)+theta,pol_rot(2,:));%將此時rot的極座標轉成xy座標(用於繪圖)
    
    [i_loc_x_top,i_loc_y_top]=pol2cart(theta+pi/2,r_rot-0.3*r);%計算頂部條棒xy座標
    [i_loc_x_bottom,i_loc_y_bottom]=pol2cart(theta-pi/2,r_rot-0.3*r);%計算底部條棒xy座標
    
    %繪出條棒電流出入紙面(有限制感應電流小於10^6視為無感應電流)
    %繪出頂部條棒
    if i_ind(i+1)>10^6
        plot(i_loc_x_top,i_loc_y_top,'ro','MarkerSize',5)
        hold on
    elseif i_ind(i+1)<-10^6
        plot(i_loc_x_top,i_loc_y_top,'rx','MarkerSize',5)
        hold on
    end
    %繪出底部條棒
    if i_ind(i+1)<-10^6
        plot(i_loc_x_bottom,i_loc_y_bottom,'co','MarkerSize',5)
        hold on
    elseif i_ind(i+1)>10^6
        plot(i_loc_x_bottom,i_loc_y_bottom,'cx','MarkerSize',5)
        hold on
    end
    
    %繪圖
    plot(0,0,'r.',xy_stat_in(1,:),xy_stat_in(2,:),'k',xy_stat_out(1,:),xy_stat_out(2,:),'k',[0,Bs_xy(1,i+1)],[0,Bs_xy(2,i+1)],'k--',xy_rot(1,:),xy_rot(2,:),'b',[0,Bw(i+1)*cos(theta)],[0,Bw(i+1)*sin(theta)])
    hold off
    axis square
    pause(0.0001);
    if i==0
        pause(5);
    end
    
end

figure(2)
subplot(3,2,1)
plot(t,ang_vel_rot)
title('rotor轉速變化圖')
xlabel('時間(s)')
ylabel('rot vel')
subplot(3,2,2)
plot(t,torque_ind)
title('感應轉矩變化圖')
xlabel('時間(s)')
ylabel('torque ind')
subplot(3,2,3)
plot(t,Bw)
title('rotor磁場變化圖')
xlabel('時間(s)')
ylabel('Bw')
subplot(3,2,4)
plot(t,e_ind)
title('感應電壓變化圖')
xlabel('時間(s)')
ylabel('e ind')
subplot(3,2,5)
plot(t,i_ind)
title('感應電流變化圖')
xlabel('時間(s)')
ylabel('i ind')
subplot(3,2,6)
plot(t,ang_accel)
title('rotor加速度變化圖')
xlabel('時間(s)')
ylabel('ang accel')