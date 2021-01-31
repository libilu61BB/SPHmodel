function [ar_x,ar_y,ae_x,ae_y] = ar_ae(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1,m_person,m_wall)
%AR_AE 计算排斥力和挤压力对各个行人产生的加速度
%   person_x 各行人粒子的x坐标
%   person_y 各行人粒子的y坐标
%   wall_x 各障碍粒子的x坐标
%   wall_y 各障碍粒子的y坐标
%   h1 常数，计算粒子密度和排斥力时使用的核半径
%   Pr_person 排斥力对各行人粒子产生的压强
%   Pr_wall 排斥力对各障碍粒子产生的压强
%   Rho_person 各行人粒子的密度
%   Rho_wall 各障碍粒子的密度outputArg1 = inputArg1;
%% 初始化参数
n = length(person_x);
s = length(wall_x);
ar_x = zeros(1,n); %ar_x 排斥力加速度在x方向上的分量
ar_y = zeros(1,n); %ar_y 排斥力加速度在y方向上的分量
ae_x = zeros(1,n); %ae_x 挤压力加速度在x方向上的分量
ae_y = zeros(1,n); %ae_y 挤压力加速度在y方向上的分量
A = 500; %排斥力压强公式参数
B = 0.001; %排斥力压强公式参数
K = 1000; %挤压力压强公式参数
h2 = 0.5; %计算障碍物排斥力所用的核半径
%% 判断坐标是否有误
if length(person_y)~=n
    error('行人的xy坐标数量不一致');
else
    if length(wall_y)~=s
        error('障碍的xy坐标数量不一致');
    end
end
%% 计算行人粒子的加速度分量
avg_Radius=mean(Radius);
Rho_p2p=m_person*(4/(pi*h1^8))*(h1^2-4*avg_Radius^2)^3+m_person*(4/(pi*h1^2)); %人与人之间的临界密度
Rho_p2w=m_person*(4/(pi*h2^8))*(h2^2-avg_Radius^2)^3+m_wall*(4/(pi*h2^2)); %人与障碍之间的临界密度
for i=1:n
    %% 计算行人与行人之间的加速度
    for j=1:n
        if j==i
            continue;
        end
        r_2 = (person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2;%计算两行人粒子之间的距离平方
        r = sqrt(r_2);
        if r_2<=h1^2 %只计算核半径范围内其他行人粒子j对粒子i密度的影响
            Rho_j2i = m_person*(4/(pi*h1^8))*(h1^2-r_2)^3+m_person*(4/(pi*h1^2)); % j对i的密度贡献加上i自身的密度
            position = (person_x(j)>=person_x(i));
            switch position
                case 1 %行人j位于行人i前方，可能产生挤压力或排斥力
                    if Rho_j2i<=Rho_p2p %未接触，计算排斥力加速度
                        Pr_j2i = A/(1+exp(-B*Rho_j2i^2));
                        abs_arj2i=m_person*(Pr_j2i/Rho_person(j)^2)*(3*(10*(h1-sqrt(r_2))^2)/(pi*h1^5));
                        ar_x(i)=ar_x(i)+abs_arj2i*(person_x(i)-person_x(j))/r; %将加速度分解到x轴
                        ar_y(i)=ar_y(i)+abs_arj2i*(person_y(i)-person_y(j))/r; %将加速度分解到y轴
                    else %接触，计算挤压力加速度
                        Pe_j2i = K*(Rho_j2i-Rho_p2p)/Rho_p2p;
                        h_ij = Radius(i)+Radius(j);
                        abs_aej2i=m_person*(Pe_j2i/Rho_person(j)^2)*(3*(10*(h_ij-sqrt(r_2))^2)/(pi*h_ij^5));
                        ae_x(i)=ae_x(i)+abs_aej2i*(person_x(i)-person_x(j))/r; %将加速度分解到x轴
                        ae_y(i)=ae_y(i)+abs_aej2i*(person_y(i)-person_y(j))/r; %将加速度分解到y轴
                    end
                case 0 %行人j位于行人i后方，只考虑挤压力
                    if Rho_j2i>Rho_p2p
                        Pe_j2i = K*(Rho_j2i-Rho_p2p)/Rho_p2p;
                        h_ij = Radius(i)+Radius(j);
                        abs_aej2i=m_person*(Pe_j2i/Rho_person(j)^2)*(3*(10*(h_ij-sqrt(r_2))^2)/(pi*h_ij^5));
                        ae_x(i)=ae_x(i)+abs_aej2i*(person_x(i)-person_x(j))/r; %将加速度分解到x轴
                        ae_y(i)=ae_y(i)+abs_aej2i*(person_y(i)-person_y(j))/r; %将加速度分解到y轴
                    end
            end
        end
    end
    %% 计算行人与障碍物之间的加速度
    d = zeros(1,s);
    for j=1:s
        d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2); %计算行人i与各个障碍物粒子之间的距离
    end
    [r,u]=min(d); %r为d中最小值，u为d中最小值的索引
    if r<=h2
        Rho_j2i = m_wall*(4/(pi*h2^8))*(h2^2-r^2)^3+m_person*(4/(pi*h2^2));
        if Rho_j2i<=Rho_p2w %人与障碍物未接触，受到排斥力
            Pr_j2i = A/(1+exp(-B*Rho_j2i^2));
            abs_arj2i=m_wall*(Pr_j2i/Rho_wall(u)^2)*(3*(10*(h2-r)^2)/(pi*h2^5));
            ar_x(i)=ar_x(i)+abs_arj2i*(person_x(i)-wall_x(u))/r; %将加速度分解到x轴
            ar_y(i)=ar_y(i)+abs_arj2i*(person_y(i)-wall_y(u))/r; %将加速度分解到y轴
        else %人与障碍物接触，受到挤压力
            Pe_j2i = K*(Rho_j2i-Rho_p2w)/Rho_p2w;
            h_ij = Radius(i);
            abs_aej2i=m_wall*(Pe_j2i/Rho_wall(u)^2)*(3*(10*(h_ij-r)^2)/(pi*h_ij^5));
            ae_x(i)=ae_x(i)+abs_aej2i*(person_x(i)-wall_x(j))/r; %将加速度分解到x轴
            ae_y(i)=ae_y(i)+abs_aej2i*(person_y(i)-wall_y(j))/r; %将加速度分解到y轴
        end
    end
end
end

