function [ar_x,ar_y,ae_x,ae_y] = ar_ae(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1,m_person,m_wall)
%AR_AE �����ų����ͼ�ѹ���Ը������˲����ļ��ٶ�
%   person_x ���������ӵ�x����
%   person_y ���������ӵ�y����
%   wall_x ���ϰ����ӵ�x����
%   wall_y ���ϰ����ӵ�y����
%   h1 ���������������ܶȺ��ų���ʱʹ�õĺ˰뾶
%   Pr_person �ų����Ը��������Ӳ�����ѹǿ
%   Pr_wall �ų����Ը��ϰ����Ӳ�����ѹǿ
%   Rho_person ���������ӵ��ܶ�
%   Rho_wall ���ϰ����ӵ��ܶ�outputArg1 = inputArg1;
%% ��ʼ������
n = length(person_x);
s = length(wall_x);
ar_x = zeros(1,n); %ar_x �ų������ٶ���x�����ϵķ���
ar_y = zeros(1,n); %ar_y �ų������ٶ���y�����ϵķ���
ae_x = zeros(1,n); %ae_x ��ѹ�����ٶ���x�����ϵķ���
ae_y = zeros(1,n); %ae_y ��ѹ�����ٶ���y�����ϵķ���
A = 500; %�ų���ѹǿ��ʽ����
B = 0.001; %�ų���ѹǿ��ʽ����
K = 1000; %��ѹ��ѹǿ��ʽ����
h2 = 0.5; %�����ϰ����ų������õĺ˰뾶
%% �ж������Ƿ�����
if length(person_y)~=n
    error('���˵�xy����������һ��');
else
    if length(wall_y)~=s
        error('�ϰ���xy����������һ��');
    end
end
%% �����������ӵļ��ٶȷ���
avg_Radius=mean(Radius);
Rho_p2p=m_person*(4/(pi*h1^8))*(h1^2-4*avg_Radius^2)^3+m_person*(4/(pi*h1^2)); %������֮����ٽ��ܶ�
Rho_p2w=m_person*(4/(pi*h2^8))*(h2^2-avg_Radius^2)^3+m_wall*(4/(pi*h2^2)); %�����ϰ�֮����ٽ��ܶ�
for i=1:n
    %% ��������������֮��ļ��ٶ�
    for j=1:n
        if j==i
            continue;
        end
        r_2 = (person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2;%��������������֮��ľ���ƽ��
        r = sqrt(r_2);
        if r_2<=h1^2 %ֻ����˰뾶��Χ��������������j������i�ܶȵ�Ӱ��
            Rho_j2i = m_person*(4/(pi*h1^8))*(h1^2-r_2)^3+m_person*(4/(pi*h1^2)); % j��i���ܶȹ��׼���i������ܶ�
            position = (person_x(j)>=person_x(i));
            switch position
                case 1 %����jλ������iǰ�������ܲ�����ѹ�����ų���
                    if Rho_j2i<=Rho_p2p %δ�Ӵ��������ų������ٶ�
                        Pr_j2i = A/(1+exp(-B*Rho_j2i^2));
                        abs_arj2i=m_person*(Pr_j2i/Rho_person(j)^2)*(3*(10*(h1-sqrt(r_2))^2)/(pi*h1^5));
                        ar_x(i)=ar_x(i)+abs_arj2i*(person_x(i)-person_x(j))/r; %�����ٶȷֽ⵽x��
                        ar_y(i)=ar_y(i)+abs_arj2i*(person_y(i)-person_y(j))/r; %�����ٶȷֽ⵽y��
                    else %�Ӵ������㼷ѹ�����ٶ�
                        Pe_j2i = K*(Rho_j2i-Rho_p2p)/Rho_p2p;
                        h_ij = Radius(i)+Radius(j);
                        abs_aej2i=m_person*(Pe_j2i/Rho_person(j)^2)*(3*(10*(h_ij-sqrt(r_2))^2)/(pi*h_ij^5));
                        ae_x(i)=ae_x(i)+abs_aej2i*(person_x(i)-person_x(j))/r; %�����ٶȷֽ⵽x��
                        ae_y(i)=ae_y(i)+abs_aej2i*(person_y(i)-person_y(j))/r; %�����ٶȷֽ⵽y��
                    end
                case 0 %����jλ������i�󷽣�ֻ���Ǽ�ѹ��
                    if Rho_j2i>Rho_p2p
                        Pe_j2i = K*(Rho_j2i-Rho_p2p)/Rho_p2p;
                        h_ij = Radius(i)+Radius(j);
                        abs_aej2i=m_person*(Pe_j2i/Rho_person(j)^2)*(3*(10*(h_ij-sqrt(r_2))^2)/(pi*h_ij^5));
                        ae_x(i)=ae_x(i)+abs_aej2i*(person_x(i)-person_x(j))/r; %�����ٶȷֽ⵽x��
                        ae_y(i)=ae_y(i)+abs_aej2i*(person_y(i)-person_y(j))/r; %�����ٶȷֽ⵽y��
                    end
            end
        end
    end
    %% �����������ϰ���֮��ļ��ٶ�
    d = zeros(1,s);
    for j=1:s
        d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2); %��������i������ϰ�������֮��ľ���
    end
    [r,u]=min(d); %rΪd����Сֵ��uΪd����Сֵ������
    if r<=h2
        Rho_j2i = m_wall*(4/(pi*h2^8))*(h2^2-r^2)^3+m_person*(4/(pi*h2^2));
        if Rho_j2i<=Rho_p2w %�����ϰ���δ�Ӵ����ܵ��ų���
            Pr_j2i = A/(1+exp(-B*Rho_j2i^2));
            abs_arj2i=m_wall*(Pr_j2i/Rho_wall(u)^2)*(3*(10*(h2-r)^2)/(pi*h2^5));
            ar_x(i)=ar_x(i)+abs_arj2i*(person_x(i)-wall_x(u))/r; %�����ٶȷֽ⵽x��
            ar_y(i)=ar_y(i)+abs_arj2i*(person_y(i)-wall_y(u))/r; %�����ٶȷֽ⵽y��
        else %�����ϰ���Ӵ����ܵ���ѹ��
            Pe_j2i = K*(Rho_j2i-Rho_p2w)/Rho_p2w;
            h_ij = Radius(i);
            abs_aej2i=m_wall*(Pe_j2i/Rho_wall(u)^2)*(3*(10*(h_ij-r)^2)/(pi*h_ij^5));
            ae_x(i)=ae_x(i)+abs_aej2i*(person_x(i)-wall_x(j))/r; %�����ٶȷֽ⵽x��
            ae_y(i)=ae_y(i)+abs_aej2i*(person_y(i)-wall_y(j))/r; %�����ٶȷֽ⵽y��
        end
    end
end
end

