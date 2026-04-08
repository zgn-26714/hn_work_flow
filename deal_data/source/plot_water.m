%% 该脚本用于绘制指定取向下水分子的gro文件
% dipole=120;%单位°
% normal=60;
clear
% OutputPath='C:\Users\Administrator\Desktop\'; % 输出地址
% filename='H2O.gro'; % 输出文件名
% twocount=1;
% twoupcount=1;
% onecount=1;
load('C:\Users\Administrator\Desktop\mycode\angle\meandata\ProDen_and_Entropy298K_2V_0ps_180normol_1dipole_ave20.mat')
i=1;j=1;
for i=1:length(angle_pro_den)
    for j=1:length(angle_pro_den{i})
len_angle=(size(angle_pro_den{i}{j},1));
H1=0;
H2=0;
for ii=1:len_angle
    for jj=1:len_angle
    pro_den=angle_pro_den{i}{j}(ii,jj)*(pi/180)^2;
%         bool=0;boll2=0;
%         if hnnormal<90 && hndipole - hnnormal > 90
%             continue
%         elseif hnnormal>90 && hndipole+hnnormal > 270
%             continue
%         end
% for dipole
%% 水分子参数
OH_length=0.09527;%单位nm
HOH_angle=104.52*pi/180;
dipolenow=ii*pi/180;
normalnow=jj*pi/180;
%% 计算原子坐标
%normal vector of H2O: L1=(x1,y1,z1); dipol vector of H2O: L2=(x2,y2,z2)
%normal vector of surface: L0=(0,0,1)
%angle between L1&L0: normal; angle between L2&L0: dipole
%norm of L1: |L1|=|OH1×OH2|=|OH1|*|OH2|*sin<OH1,OH2> 
%|OH|=0.9527,<OH1,OH2>=104.52°
norm_L1=power(OH_length,2)*sin(HOH_angle);
%cos(normal)=(L1·L0)/(|L1|*|L0|),and assume x1=0, then
x1=0;
y1=norm_L1*sin(normal);
z1=norm_L1*cos(normal);
% now we get L1=(0,y1,z1). As L1 is the normal vector of H2O and
% L2 is the dipol vector of H2O, L1 is normal to L2
% so y1y2+z1z2=0   [1]
% norm of L2: |L2|=|OH1+OH2|=sqrt(|OH1|^2+|OH2|^2+2|a|*|b|*cos<OH1,OH2>)
% while |L2|=sqrt(x2^2+y2^2+z2^2)   [2]
norm_L2=sqrt(power(OH_length,2)+power(OH_length,2)+2*power(OH_length,2)*cos(HOH_angle));
% cos(dipole)=(L2·L0)/(|L2|*|L0|), i.e. cos(dipole)=z2/|L2|, then
z2=norm_L2*cos(dipole);% [3]
% simultaneous formulas [1] [2] [3], we get
y2=-z1*z2/y1;
x2=sqrt(power(norm_L2,2)-power(y2,2)-power(z2,2));
L1=[x1,y1,z1];L2=[x2,y2,z2];
% now we get L2=(x2,y2,z2)
% assume OH vecotr = (a,b,c),then
% cos(HOH_angle*0.5)=OH·L2/|OH||L2|    [4]
% OH·L1=y1*b+z1*c=0    [5]
% |OH|=sqrt(a^2+b^2+c^2)    [6]
% simultaneous formulas [4] [5] [6], we get
syms B;
Coe1=power(((y1*z2-y2*z1)/(z1*x2)),2)+1+power(y1/z1,2);
Coe2=2*((y1*z2-y2*z1)/(z1*x2))*OH_length*norm_L2*cos(HOH_angle*0.5)/x2;
Coe3=power(OH_length*norm_L2*cos(HOH_angle*0.5)/x2,2)-power(OH_length,2); 
Formula=B^2*Coe1+B*Coe2+Coe3;
B=solve(Formula==0);
b=double(B);
c=-b*y1/z1;
a=(OH_length*norm_L2*cos(HOH_angle*0.5)-b*y2-c*z2)/x2;
norm_OH=zeros(1,2);
try
    for i=[1,2]
        norm_OH(i)=sqrt(a(i)^2+b(i)^2+c(i)^2);
        a(i)=a(i)*OH_length/norm_OH(i);
        b(i)=b(i)*OH_length/norm_OH(i);
        c(i)=c(i)*OH_length/norm_OH(i);
    end
catch 
    continue
end

% assume coordinate of OW: (1,1,1)
% coordinate of HW1: (a(1),b(1),c(1)+1.5);coordinate of HW2: (a(2),b(2),c(2)+1.5);
X=[a(1)+0.25,0.25,a(2)+0.25];Y=[b(1)+0.25,0.25,b(2)+0.25];Z=[c(1)+0.25,0.25,c(2)+0.25];
H1=Z(1)+H1*pro_den;H2=Z(3)+H2*pro_den;
% vector1=O-H1;
% vector2=O-H2;
% bool1=sum(vector1.*[0 0 1]) < 0 && sum(vector2.*[0 0 1]) < 0;
% bool2=sum(vector1.*[0 0 1]) > 0 && sum(vector2.*[0 0 1]) > 0;
% if bool1
%     twodownnormal(twocount)=hnnormal;
%     twodowndipole(twocount)=hndipole;
%     twocount=twocount+1;
% elseif bool2
%     twoupnormal(twoupcount)=hnnormal;
%     twoupdipole(twoupcount)=hndipole;
%     twoupcount=twoupcount+1;
% else
%     onedownnormal(onecount)=hnnormal;
%     onedowndipole(onecount)=hndipole;
%     onecount=onecount+1;
% end
    end
end
H1all(i,j)=H1/len_angle^2;
H2all(i,j)=H2/len_angle^2;
    end
end
% for hndipole=dipole{i}{j}
%     for hnnormal=normal{i}{j}
%         a=hndipole;
%         b=hnnormal;
%     end
% end
