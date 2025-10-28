global RBmap


RBmap = zeros(101,3);
RBmap(51,:) = 1;

RBmap(1:50,:) = [linspace(0,1,50)' linspace(0,1,50)' ones(50,1)];
RBmap(52:101,:) = [ones(50,1) linspace(1,0,50)' linspace(1,0,50)'];

function plot_surf(xindex,yindex,zdata,vi)
    global RBmap
    surf(xindex,yindex,zdata)
    hold on;
    % plot3([1 1],[0 20],[10 10],'k--','LineWidth',1)
    % plot3([0.5 0.5],[0 20],[10 10],'k:','LineWidth',1)
    shading interp
    view(vi)
    colormap(RBmap)
    % xticks([0 1 2])
    % xticklabels([])
    % ylim([0 30])
    % yticks([0 15 30])
    % yticklabels([0 15 30])
    % caxis([-0.04 0.04])
    % ylabel('{\itt } (ps)')
    % xlim([0 2])
    % title('negative')
    bar=colorbar();
    set(bar,'Units','centimeters','Position',[14 4 0.5 8])
    set(gca,'Units','centimeters','Position',[2 4 10 8])
end

function remake_angle(data,angle1Bin,angle2Bin,nz)
    % angle1Bin;%normal
    % angle2Bin;%dipole
    for i = 1 : size(data,1)/nz
        tmp = data( (i-1)*nz+1 : i*nz, :);
        for j = 1:nz
            angle_prby{i}{j} = reshape(tmp(j,:), [angle2Bin, angle1Bin]);
            if(sum(sum(angle_prby{i}{j}))~=0)
                angle_prby{i}{j} = angle_prby{i}{j}/sum(sum(angle_prby{i}{j}));
            end
        end
    end
    lz =10;
    
end

function result = to_average(angle_prby,layer,lz)
    tmp3D = to3D(angle_prby,layer,lz);
    % tmp3D = fliplr(tmp3D);
    for i = 1: size(tmp3D,1)     % 时间维度
        result(i) = 0;
        for j = 1:size(tmp3D,2)
            result(i) = tmp3D(i,j) * (j-0.5)/size(tmp3D,2) * 180 + result(i);
        end
    end
end
    
    
function neg3D=to3D(angle_prby,layer,lz)
    negAng = divide_layer(angle_prby,layer,lz);
    if(~ismember(1,size(negAng{1})))
        disp("there is error data!");
        exit
    end
    for i = 1:length(negAng)
        neg3D(i,:) = negAng{i};
    end
end
    
function layer_angle = divide_layer(ori_data, layer, lz)
    dz = lz / length(ori_data{1});
    begin = int32(layer(1) / dz);
    if (0 == begin)
        begin = 1;
    end
    last = int32 (layer(2) / dz);
    for t = 1:length(ori_data)
        tmp{t} = zeros(size(ori_data{t}{1}));
        for i = begin : last
            tmp{t} = tmp{t} + ori_data{t}{i};
        end
        tmp{t} = tmp{t} / double(last - begin + 1);
    end
    layer_angle = tmp;
end




