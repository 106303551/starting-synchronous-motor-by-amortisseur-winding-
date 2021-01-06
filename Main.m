clc
clear

%�q���ѼƳ]�w
r_stat_in=15;%stator����b�|
r_stat_out=17;%stator�~��b�|
width_rot=2;%rotor�ӳ��e��
length_rot=14;%rotor�ӳ�����
r=1.5;%rotor�b��b�|
x=r*cosd(0:1:180);%rotor�b�ꪺx�b���q
y=r*sind(0:1:180);%rotor�b�ꪺy�b���q
r_rot=length_rot/2+r;%rotor�b�|
theta=0;%rotor��l����

%stator�]�w
Bs=12;%stator�ϳ��¶q
ang_freq=2*pi;%���W�v(�]�w�w���ϳ���@��@��)
dt=0.01%�ɶ����q����
t=0:dt:25;%�`�ɶ�����

%�]�wstator�y��
xy_stat_in=[r_stat_in*cosd(0:4:360);r_stat_in*sind(0:4:360)];%stator����xy�y��
xy_stat_out=[r_stat_out*cosd(0:4:360);r_stat_out*sind(0:4:360)];%stator�~��xy�y��

%�]�wrotor��xy�y��
xy_rot=[[width_rot/2,width_rot/2,width_rot/2+0.5,x,-width_rot/2,-width_rot/2,-width_rot/2,-x,width_rot/2,width_rot/2];
    [0,length_rot/2,length_rot/2,length_rot/2+y,length_rot/2,0,-length_rot/2,-length_rot/2-y,-length_rot/2,0]];
[pol_rot(1,:),pol_rot(2,:)]=cart2pol(xy_rot(1,:),xy_rot(2,:));%�Nrotor��xy�b�y���ন���y��(�Ω��K�p��)
Bs_x=Bs*cos(ang_freq*t+pi/2);%stator�ϳ���x���q
Bs_y=Bs*sin(ang_freq*t+pi/2);%stator�ϳ���y���q
Bs_xy=[Bs_x;Bs_y];%stator�ϳ�xy���q�x�}


%�ѼƳ]�w
R_winding=2*10^(-3);%�����u¶�չq��
G=1;%�X��Y��
mu=1.2566*10^(-6);%�ɺϲv(��)
k=r_rot*r_rot*pi*G/mu;%�q�����c�]�� k=A(r^2*pi)*G/mu
moment_inertia=10^10;%��ʺD�q

%�]�w�Ω��x�s�H�ɶ����ܪ����z�q���x�}

e_ind=zeros(1,size(t,2));%�P���q��
i_ind=zeros(1,size(t,2));%�P���q�y
Bw=zeros(1,size(t,2));%�ೡ�ϳ�
torque_ind=zeros(1,size(t,2));%�P����x
ang_accel=zeros(1,size(t,2));%���t��
ang_vel_rot=zeros(1,size(t,2));%rotor��t

%�p��|�H�ɶ����ܪ����z�q(�Hdt=0.01���p�⪺�ɶ����j)
for i=0:1:size(t,2)-1 %dt������
    %i+1���ثe��U�����z�q i+2���U��dt�����z�q
    [vel_x,vel_y]=pol2cart(theta,r_rot*(ang_freq-ang_vel_rot(i+1)));%�p����ά۹�ϳ��bxandy�W���t�פ��q
    e_ind(i+1)=(vel_x*Bs_xy(2,i+1)-vel_y*Bs_xy(1,i+1))*length_rot;%�P���q��=(V�~�nB)*l(1*2�x�}�M1*2�x�}�~�n�Y����e�ۭ�)
    i_ind(i+1)=e_ind(i+1)/R_winding;%�p��P���q�y
    Bw(i+1)=mu*i_ind(i+1)/G;%�p��¶�պϳ����¶q
    torque_ind(i+1)=k*((Bw(i+1)*cos(theta))*Bs_xy(2,i+1)-(Bw(i+1)*sin(theta))*Bs_xy(1,i+1));%�P����x=k*(¶�պϳ��~�nstat�ϳ�)(��e�ۭ�)
    ang_accel(i+1)=torque_ind(i+1)/moment_inertia;%���[�t��
    
    if i<size(t,2)-1 %���קK�h�p��@��dt����t �Φ�������(�T�Ot and ang_vel_rot�x�}���פ@��)
        ang_vel_rot(i+2)=ang_vel_rot(i+1)+ang_accel(i+1)*dt;%�p���H��dt���ܪ���t(�U��dt����t����l��t+(���[�t��*dt))
    end
    
    theta_plus=ang_vel_rot(i+1)*dt;%�p��o��dt���ܪ�����(���b�o��dt������t���`��)
    theta=theta+theta_plus;%�p��U��dt�ɪ�����(��dt�ɪ�����+��dt���ܪ�����)
    [xy_rot(1,:),xy_rot(2,:)]=pol2cart(pol_rot(1,:)+theta,pol_rot(2,:));%�N����rot�����y���নxy�y��(�Ω�ø��)
    
    [i_loc_x_top,i_loc_y_top]=pol2cart(theta+pi/2,r_rot-0.3*r);%�p�⳻������xy�y��
    [i_loc_x_bottom,i_loc_y_bottom]=pol2cart(theta-pi/2,r_rot-0.3*r);%�p�⩳������xy�y��
    
    %ø�X���ιq�y�X�J�ȭ�(������P���q�y�p��10^6�����L�P���q�y)
    %ø�X��������
    if i_ind(i+1)>10^6
        plot(i_loc_x_top,i_loc_y_top,'ro','MarkerSize',5)
        hold on
    elseif i_ind(i+1)<-10^6
        plot(i_loc_x_top,i_loc_y_top,'rx','MarkerSize',5)
        hold on
    end
    %ø�X��������
    if i_ind(i+1)<-10^6
        plot(i_loc_x_bottom,i_loc_y_bottom,'co','MarkerSize',5)
        hold on
    elseif i_ind(i+1)>10^6
        plot(i_loc_x_bottom,i_loc_y_bottom,'cx','MarkerSize',5)
        hold on
    end
    
    %ø��
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
title('rotor��t�ܤƹ�')
xlabel('�ɶ�(s)')
ylabel('rot vel')
subplot(3,2,2)
plot(t,torque_ind)
title('�P����x�ܤƹ�')
xlabel('�ɶ�(s)')
ylabel('torque ind')
subplot(3,2,3)
plot(t,Bw)
title('rotor�ϳ��ܤƹ�')
xlabel('�ɶ�(s)')
ylabel('Bw')
subplot(3,2,4)
plot(t,e_ind)
title('�P���q���ܤƹ�')
xlabel('�ɶ�(s)')
ylabel('e ind')
subplot(3,2,5)
plot(t,i_ind)
title('�P���q�y�ܤƹ�')
xlabel('�ɶ�(s)')
ylabel('i ind')
subplot(3,2,6)
plot(t,ang_accel)
title('rotor�[�t���ܤƹ�')
xlabel('�ɶ�(s)')
ylabel('ang accel')