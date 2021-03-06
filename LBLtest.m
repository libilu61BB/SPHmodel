clear;
%% 设置障碍物坐标、行人坐标和出口坐标
% 原设置orz
% 15m×15m正方形空间，出口宽度2m
wall_x1=(15:-0.1:0);wall_y1=zeros(1,length(wall_x1));
wall_y2=(0:0.1:15);wall_x2=zeros(1,length(wall_y2));
wall_x3=(0:0.1:15);wall_y3=15*ones(1,length(wall_x3));
wall_y4=(15:-0.1:8.5);wall_x4=15*ones(1,length(wall_y4));
wall_y5=(6.5:-0.1:0);wall_x5=15*ones(1,length(wall_y5));
wall_x=[wall_x5 wall_x1 wall_x2 wall_x3 wall_x4];
wall_y=[wall_y5 wall_y1 wall_y2 wall_y3 wall_y4];
person_x = 9*rand(1,100)+1;
person_y = 9*rand(1,100)+1;
exit_x = 20;
exit_y = 7.5;
end_x = 16;

% 啊 让我们来康康2m×100m的走廊呢~
% wall_x1 = (-10:0.1:100);
% wall_y1 = zeros(1, length(wall_x1));
% wall_x2 = (-100:0.1:100);
% wall_y2 = 2 * ones(1, length(wall_x2));
% wall_x = [wall_x1 wall_x2];
% wall_y = [wall_y1 wall_y2];
% %在空间内随机生成点用于模拟行人
% % person_x=zeros(1,100)+10;
% % person_y=2*rand(1,100);
% % person_x=[0.3:0.3:30];
% % person_y=1.4*rand(1,100)+0.3;
% load personIni.mat
% exit_x=1000;%出口x坐标
% exit_y=1;%出口y坐标
% end_x = 100;%行人消失的点

%% 计算坐标，绘制图像
n=length(person_x);
s=length(wall_x);
h1=5;%计算密度和排斥力时使用的核半径
Radius=0.25*ones(1,n);%假设行人的半径均为0.25m
m_person=70;%行人的质量
m_wall=500;%障碍物的质量
a_max = 5; %行人的主动力加速度上限
v0=2; %行人的期望速度值
u=2; %粘度，用于计算粒子之间摩擦力产生的加速度
vx=zeros(1,n);%行人速度在x方向上的分量，初始时刻为0
vy=zeros(1,n);%行人速度在y方向上的分量，初始时刻为0
t1 = 0.25;%人的反应时间，s
T=200; %模拟总时间
sum_escape=0;%统计已疏散的人数
P=1;%熟悉逃生路线的行人比例
P_f=0.8;%从众程度
dt=0.03;

for t=0:dt:T
    %% 计算排斥力和挤压力产生的加速度
    [Rho_person,Rho_wall]=density(person_x,person_y,wall_x,wall_y,h1);%调用函数density计算t时刻的密度
%     [ar_x,ar_y]=a_repul(person_x,person_y,wall_x,wall_y,Rho_person ,Rho_wall,Radius,h1);%调用函数a_repul计算排斥力产生的加速度
%     [ae_x,ae_y]=a_extru(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1);%调用函数a_extru计算挤压力产生的加速度
    [ar_x,ar_y,ae_x,ae_y] = ar_ae(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1,m_person,m_wall); %计算加速度
    
    %% 计算摩擦力产生的加速度
    av_x=zeros(1,n);%初始化摩擦力产生的x方向加速度
    av_y=zeros(1,n);%初始化摩擦力产生的y方向加速度
    for i=1:n
        for j=1:n %计算行人与行人之间的摩擦力产生的加速度
            if j==i
                continue;
            end
            r=sqrt((person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2); %两行人粒子之间的距离
            if r<=(Radius(i)+Radius(j))
                av_x(i)=av_x(i)+u*m_person*((vx(j)-vx(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
                av_y(i)=av_y(i)+u*m_person*((vy(j)-vy(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
            end
        end
        d=zeros(1,s);%初始化行人粒子i与各障碍粒子的距离
        for j=1:s
            d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
        end
        [r,m]=min(d); %r为d中最小值，m为d中最小值的索引
        if r<=Radius(i) %行人与障碍之间的摩擦力产生的加速度
            av_x(i)=av_x(i)+u*m_wall*((0-vx(i))/(Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
            av_y(i)=av_y(i)+u*m_wall*((0-vy(i))/(Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
        end
    end
    
    %% 计算行人主动力产生的加速度
    e_x=zeros(1,n); %初始化方向向量的x坐标
    e_y=zeros(1,n); %初始化方向向量的y坐标
    am_x=zeros(1,n); %初始化主动力产生的x方向加速度
    am_y=zeros(1,n); %初始化主动力产生的y方向加速度
    for i=1:(P*n) %计算熟悉逃生路线的行人的运动方向
        r=sqrt((person_x(i)-exit_x)^2+(person_y(i)-exit_y)^2); %熟悉逃生路线的行人粒子与出口之间的距离
        e_x(i)=(exit_x-person_x(i))/r;
        e_y(i)=(exit_y-person_y(i))/r;
    end
    for i=(P*n+1):n %为不熟悉逃生路线的行人产生一个随机的运动方向
        e_x(i)=-1+2*rand;
        e_y(i)=-1+2*rand;
        r=sqrt(e_x(i)^2+e_y(i)^2);%方向向量的模
        e_x(i)=e_x(i)/r;%将方向向量转换为单位方向向量
        e_y(i)=e_y(i)/r;%将方向向量转换为单位方向向量
    end
    for i=(P*n+1):n %计算不熟悉逃生路线的行人在从众行为影响下的运动方向
        e_x(i)=(1-P_f)*e_x(i);
        e_y(i)=(1-P_f)*e_y(i);
        for j=1:n
            if j==i
                continue;
            end
            r_2=(person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2; %两行人粒子距离的平方
            if r_2<=h1^2
                e_x(i)=e_x(i)+P_f*(m_person*e_x(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r_2)^3;
                e_y(i)=e_y(i)+P_f*(m_person*e_y(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r_2)^3;
            end
        end
        r=sqrt(e_x(i)^2+e_y(i)^2);%方向向量的模
        e_x(i)=e_x(i)/r;%将方向向量转换为单位方向向量
        e_y(i)=e_y(i)/r;%将方向向量转换为单位方向向量
        r=sqrt((person_x(i)-exit_x)^2+(person_y(i)-exit_y)^2);
        if r<=5 %当不熟悉疏散路线的行人接近出口时，将其视为熟悉疏散路线的行人
            if person_x(i)<=exit_x
                e_x(i)=(exit_x-person_x(i))/r;
                e_y(i)=(exit_y-person_y(i))/r;
            end
        end
    end
    for i=1:n
%       am_x(i)=(v0*e_x(i)-vx(i))/dt;
%       am_y(i)=(v0*e_y(i)-vy(i))/dt;
        am_x(i)=(v0*e_x(i)-vx(i))/t1;%设置主动力加速度的上限
        am_y(i)=(v0*e_y(i)-vy(i))/t1;
    end
    
    %% 计算行人的位置
    ax=am_x+ar_x+ae_x+av_x;%1行n列，t时刻各行人粒子x方向的合加速度
    ay=am_y+ar_y+ae_y+av_y;%1行n列，t时刻各行人粒子y方向的合加速度
%     for i=1:n
%         ar = sqrt(ax(i)^2+ay(i)^2);
%         if ar>a_max
%             ax(i) = ax(i)*a_max/ar;
%             ay(i) = ay(i)*a_max/ar;
%         end
%     end   
    vx=vx+ax*dt; %计算下一时刻的x方向速度
    vy=vy+ay*dt; %计算下一时刻的y方向速度 
    for i=1:n %若下一时刻的速度大于v0，则将其缩小到v0
        vr=sqrt(vx(i)^2+vy(i)^2);
        if vr>v0
            vx(i)=vx(i)*v0/vr;
            vy(i)=vy(i)*v0/vr;
        end
    end
    person_x=person_x+vx*dt; %使用dt时间段内的平均速度计算x方向的位移
    person_y=person_y+vy*dt; %使用dt时间段内的平均速度计算y方向的位移

    %% 绘制图像
    for i=1:n %统计逃离人数
        if person_x(i)>(end_x)
            sum_escape=sum_escape+1;
            person_x(i)=nan;
            person_y(i)=nan;
        end
    end
    
%     plot(wall_x1,wall_y1,'LineWidth',2,'Color','k');
%     hold on;
%     plot(wall_x2,wall_y2,'LineWidth',2,'Color','k');
%     plot(person_x,person_y,'.','MarkerSize',15)
%     axis([-1 100 -1 3]);%设置显示范围:2×100m通道
    
    plot(wall_x,wall_y)
    hold on
    plot(person_x,person_y,'.','MarkerSize',15)
    axis([-1 16 -1 16]);%设置显示范围:15×15m房间
    
    str_time=sprintf('疏散时间：%.2f',t);
    str_escape=sprintf('疏散人数：%.0f',sum_escape);
    text(max(xlim)*0.5-10,-0.5,str_time);
    text(max(xlim)*0.5-10,-0.7,str_escape);
    axis on;
    hold off;
    pause(0.001);
    if sum_escape>=n
        break;
    end
    
    %% 制作GIF
    %     frame = getframe(gcf);
    %     im = frame2im(frame);
    %     [I,map] = rgb2ind(im,256);
    %     if t == 0
    %         imwrite(I,map,'疏散模拟.gif','gif','Loopcount',inf,'DelayTime',0.01);
    %     else
    %         imwrite(I,map,'疏散模拟.gif','gif','WriteMode','append','DelayTime',0.01);
    %     end

end            