%% 该脚本用于判断两个H在O的哪一侧
clear

% circle=0.3;c_n=1;diedai=0.3;
count=ones(180)*-1;
for i=1:180
    i
    for j=1:180
        inormal=i;
        idipole=j;
        filename=sprintf('%g_%g.gro',inormal,idipole);
        
        %% 水分子参数
        OH_length=0.09527;%单位nm
        HOH_angle=104.52*pi/180;
        dipole=idipole*pi/180;
        normal=inormal*pi/180;
        %% 计算原子坐标
        norm_L1=power(OH_length,2)*sin(HOH_angle);

        x1=0;
        y1=norm_L1*sin(normal);
        z1=norm_L1*cos(normal);
 
        norm_L2=sqrt(power(OH_length,2)+power(OH_length,2)+2*power(OH_length,2)*cos(HOH_angle));

        z2=norm_L2*cos(dipole);% [3]

        y2=-z1*z2/y1;
        x2=sqrt(power(norm_L2,2)-power(y2,2)-power(z2,2));
        L1=[x1,y1,z1];L2=[x2,y2,z2];

        syms B;
        Coe1=power(((y1*z2-y2*z1)/(z1*x2)),2)+1+power(y1/z1,2);
        Coe2=2*((y1*z2-y2*z1)/(z1*x2))*OH_length*norm_L2*cos(HOH_angle*0.5)/x2;
        Coe3=power(OH_length*norm_L2*cos(HOH_angle*0.5)/x2,2)-power(OH_length,2); 
        Formula=B^2*Coe1+B*Coe2+Coe3;
        lastwarn('');
        B=solve(Formula==0);
        warningMessage = lastwarn;
        if ~isempty(warningMessage)
            continue;
        end
        b=double(B);
        c=-b*y1/z1;
        a=(OH_length*norm_L2*cos(HOH_angle*0.5)-b*y2-c*z2)/x2;
        norm_OH=zeros(1,2);
        for ii=[1,2]
            norm_OH(ii)=sqrt(a(ii)^2+b(ii)^2+c(ii)^2);
            a(ii)=a(ii)*OH_length/norm_OH(ii);
            b(ii)=b(ii)*OH_length/norm_OH(ii);
            c(ii)=c(ii)*OH_length/norm_OH(ii);
        end
        if(~isreal(a(1)) || ~isreal(a(2)))
            continue;
        end
        if(~isreal(b(1)) || ~isreal(b(2)))
            continue;
        end
        if(~isreal(c(1)) || ~isreal(c(2)))
            continue;
        end
        
        if ( (c(1) <= 0 && c(2) <=0) || (c(1) >= 0 && c(2) >=0))
            count(i,j)=1;
        else
            count(i,j)=0;
        end
    end
end
