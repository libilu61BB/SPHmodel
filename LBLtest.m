clear;
%% �����ϰ������ꡢ��������ͳ�������
% ԭ����orz
%15m��15m�����οռ䣬���ڿ��2m
% wall_x1=(15:-0.1:0);wall_y1=zeros(1,length(wall_x1));
% wall_y2=(0:0.1:15);wall_x2=zeros(1,length(wall_y2));
% wall_x3=(0:0.1:15);wall_y3=15*ones(1,length(wall_x3));
% wall_y4=(15:-0.1:8.5);wall_x4=15*ones(1,length(wall_y4));
% wall_y5=(6.5:-0.1:0);wall_x5=15*ones(1,length(wall_y5));
% wall_x=[wall_x5 wall_x1 wall_x2 wall_x3 wall_x4];
% wall_y=[wall_y5 wall_y1 wall_y2 wall_y3 wall_y4];
% person_x = 9*rand(1,100)+1;
% person_y = 9*rand(1,100)+1;
% exit_x = 20;
% exit_y = 7.5;
% end_x = 13;

% % �� ������������2m��100m��������~
wall_x1 = (-100:0.1:100);
wall_y1 = zeros(1, length(wall_x1));
wall_x2 = (-100:0.1:100);
wall_y2 = 2 * ones(1, length(wall_x2));
wall_x = [wall_x1 wall_x2];
wall_y = [wall_y1 wall_y2];
%�ڿռ���������ɵ�����ģ������
% person_x=zeros(1,100)+10;
% person_y=2*rand(1,100);
person_x=rand(1,100)+30;
person_y=1.5*rand(1,100)+0.25;
exit_x=1000;%����x����
exit_y=1;%����y����
end_x = 100;%������ʧ�ĵ�

%% �������꣬����ͼ��
n=length(person_x);
s=length(wall_x);
h1=5;%�����ܶȺ��ų���ʱʹ�õĺ˰뾶
Radius=0.25*ones(1,n);%�������˵İ뾶��Ϊ0.25m
m_person=70;%���˵�����
m_wall=500;%�ϰ��������
a_max = 3; %���˵ļ��ٶ�����
v0=2; %���˵������ٶ�ֵ
u=2; %ճ�ȣ����ڼ�������֮��Ħ���������ļ��ٶ�
vx=zeros(1,n);%�����ٶ���x�����ϵķ�������ʼʱ��Ϊ0
vy=zeros(1,n);%�����ٶ���y�����ϵķ�������ʼʱ��Ϊ0
T=200; %ģ����ʱ��
sum_escape=0;%ͳ������ɢ������
P=1;%��Ϥ����·�ߵ����˱���
P_f=0.8;%���ڳ̶�
dt=0.03;

for t=0:dt:T
    %% �����ų����ͼ�ѹ�������ļ��ٶ�
    [Rho_person,Rho_wall]=density(person_x,person_y,wall_x,wall_y,h1);%���ú���density����tʱ�̵��ܶ�
%     [ar_x,ar_y]=a_repul(person_x,person_y,wall_x,wall_y,Rho_person ,Rho_wall,Radius,h1);%���ú���a_repul�����ų��������ļ��ٶ�
%     [ae_x,ae_y]=a_extru(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1);%���ú���a_extru���㼷ѹ�������ļ��ٶ�
    [ar_x,ar_y,ae_x,ae_y] = ar_ae(person_x,person_y,wall_x,wall_y,Rho_person,Rho_wall,Radius,h1,m_person,m_wall); %������ٶ�
    
    %% ����Ħ���������ļ��ٶ�
    av_x=zeros(1,n);%��ʼ��Ħ����������x������ٶ�
    av_y=zeros(1,n);%��ʼ��Ħ����������y������ٶ�
    for i=1:n
        for j=1:n %��������������֮���Ħ���������ļ��ٶ�
            if j==i
                continue;
            end
            r=sqrt((person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2); %����������֮��ľ���
            if r<=(Radius(i)+Radius(j))
                av_x(i)=av_x(i)+u*m_person*((vx(j)-vx(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
                av_y(i)=av_y(i)+u*m_person*((vy(j)-vy(i))/(Rho_person(i)*Rho_person(j)))*(1800*(Radius(i)+Radius(j)-r))/(45*pi*(Radius(i)+Radius(j))^5);
            end
        end
        d=zeros(1,s);%��ʼ����������i����ϰ����ӵľ���
        for j=1:s
            d(j)=sqrt((person_x(i)-wall_x(j))^2+(person_y(i)-wall_y(j))^2);
        end
        [r,m]=min(d); %rΪd����Сֵ��mΪd����Сֵ������
        if r<=Radius(i) %�������ϰ�֮���Ħ���������ļ��ٶ�
            av_x(i)=av_x(i)+u*m_wall*((0-vx(i))/(Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
            av_y(i)=av_y(i)+u*m_wall*((0-vy(i))/(Rho_person(i)*Rho_wall(m)))*(1800*(Radius(i)-r))/(45*pi*(Radius(i))^5);
        end
    end
    
    %% �������������������ļ��ٶ�
    e_x=zeros(1,n); %��ʼ������������x����
    e_y=zeros(1,n); %��ʼ������������y����
    am_x=zeros(1,n); %��ʼ��������������x������ٶ�
    am_y=zeros(1,n); %��ʼ��������������y������ٶ�
    for i=1:(P*n) %������Ϥ����·�ߵ����˵��˶�����
        r=sqrt((person_x(i)-exit_x)^2+(person_y(i)-exit_y)^2); %��Ϥ����·�ߵ��������������֮��ľ���
        e_x(i)=(exit_x-person_x(i))/r;
        e_y(i)=(exit_y-person_y(i))/r;
    end
    for i=(P*n+1):n %Ϊ����Ϥ����·�ߵ����˲���һ��������˶�����
        e_x(i)=-1+2*rand;
        e_y(i)=-1+2*rand;
        r=sqrt(e_x(i)^2+e_y(i)^2);%����������ģ
        e_x(i)=e_x(i)/r;%����������ת��Ϊ��λ��������
        e_y(i)=e_y(i)/r;%����������ת��Ϊ��λ��������
    end
    for i=(P*n+1):n %���㲻��Ϥ����·�ߵ������ڴ�����ΪӰ���µ��˶�����
        e_x(i)=(1-P_f)*e_x(i);
        e_y(i)=(1-P_f)*e_y(i);
        for j=1:n
            if j==i
                continue;
            end
            r_2=(person_x(i)-person_x(j))^2+(person_y(i)-person_y(j))^2; %���������Ӿ����ƽ��
            if r_2<=h1^2
                e_x(i)=e_x(i)+P_f*(m_person*e_x(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r_2)^3;
                e_y(i)=e_y(i)+P_f*(m_person*e_y(j)/Rho_person(j))*(4/(pi*h1^8))*(h1^2-r_2)^3;
            end
        end
        r=sqrt(e_x(i)^2+e_y(i)^2);%����������ģ
        e_x(i)=e_x(i)/r;%����������ת��Ϊ��λ��������
        e_y(i)=e_y(i)/r;%����������ת��Ϊ��λ��������
        r=sqrt((person_x(i)-exit_x)^2+(person_y(i)-exit_y)^2);
        if r<=5 %������Ϥ��ɢ·�ߵ����˽ӽ�����ʱ��������Ϊ��Ϥ��ɢ·�ߵ�����
            if person_x(i)<=exit_x
                e_x(i)=(exit_x-person_x(i))/r;
                e_y(i)=(exit_y-person_y(i))/r;
            end
        end
    end
    for i=1:n
%       am_x(i)=(v0*e_x(i)-vx(i))/dt;
%       am_y(i)=(v0*e_y(i)-vy(i))/dt;
        am_x(i)=min((v0*e_x(i)-vx(i))/dt,a_max);%�������������ٶȵ�����
        am_y(i)=min((v0*e_y(i)-vy(i))/dt,a_max);
    end
    
    %% �������˵�λ��
    ax=min(a_max,am_x+ar_x+ae_x+av_x);%1��n�У�tʱ�̸���������x����ĺϼ��ٶ�
    ay=min(a_max,am_y+ar_y+ae_y+av_y);%1��n�У�tʱ�̸���������y����ĺϼ��ٶ�
    old_vx=vx;%����dtʱ���ڵ�x����ƽ���ٶ�
    old_vy=vy;%����dtʱ���ڵ�y����ƽ���ٶ�
    vx=vx+ax*dt; %������һʱ�̵�x�����ٶ�
    vy=vy+ay*dt; %������һʱ�̵�y�����ٶ�
    for i=1:n %����һʱ�̵��ٶȴ���v0��������С��v0
        vr=sqrt(vx(i)^2+vy(i)^2);
        if vr>=v0
            vx(i)=vx(i)*v0/vr;
            vy(i)=vy(i)*v0/vr;
        end
    end
%     for i=1:n
%         if person_x(i)>100 %����ͨ�����ں���x�����뿪
%             vx(i)=v0;
%             vy(i)=0;
%         end
%     end
    avg_vx=(old_vx+vx)/2;
    avg_vy=(old_vy+vy)/2; 
    person_x=person_x+avg_vx*dt; %ʹ��dtʱ����ڵ�ƽ���ٶȼ���x�����λ��
    person_y=person_y+avg_vy*dt; %ʹ��dtʱ����ڵ�ƽ���ٶȼ���y�����λ��

    %% ����ͼ��
    for i=1:n %ͳ����������
        if person_x(i)>(end_x)
            sum_escape=sum_escape+1;
            person_x(i)=nan;
            person_y(i)=nan;
        end
    end
    plot(wall_x1,wall_y1,'LineWidth',2,'Color','k');
    hold on;
    plot(wall_x2,wall_y2,'LineWidth',2,'Color','k');
%     plot(wall_x,wall_y)
%     hold on
    plot(person_x,person_y,'.')
    axis([-1 100 -1 3]);%������ʾ��Χ
    str_time=sprintf('��ɢʱ�䣺%.2f',t);
    str_escape=sprintf('��ɢ������%.0f',sum_escape);
    text(max(xlim)*0.5-10,-0.5,str_time);
    text(max(xlim)*0.5-10,-0.7,str_escape);
    axis on;
    hold off;
    pause(0.001);
    if sum_escape>=n
        break;
    end
    
    %% ����GIF
    %     frame = getframe(gcf);
    %     im = frame2im(frame);
    %     [I,map] = rgb2ind(im,256);
    %     if t == 0
    %         imwrite(I,map,'��ɢģ��.gif','gif','Loopcount',inf,'DelayTime',0.01);
    %     else
    %         imwrite(I,map,'��ɢģ��.gif','gif','WriteMode','append','DelayTime',0.01);
    %     end

end            