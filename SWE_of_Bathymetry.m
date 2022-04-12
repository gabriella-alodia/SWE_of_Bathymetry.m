%% Slope-weighted eccentricity (SWE)
%  Alodia et al. (2022)
%  A geomorphometric algorithm to obtain the numerical description of both 
%  magmatic and tectonic crust in a slow-spreading ridge through a series 
%  of calculation based on the distribution of the azimuth and plunge 
%  observed in the seafloor morphology.

%% Input needed
%  A gridded shipborne multibeam bathymetry (depths in metres)
%  in *.xyz format (here: 'Input_Bathymetry_15s.xyz')

%% Output
%  1. Terrain eccentricity (here: 'Output_eccentricity.xyz')
%  2. Weight matrix: 1-sin(slope) (here: 'Output_weight.xyz')
%  3. SWE: Slope-weighted eccentricity (here: 'Output_SWE.xyz')
%  4. Masked SWE (here: 'Output_SWE_masked.xyz')
%  Each output is exported in *.xyz format. An explanation on converting
%  *.xyz files into *.grd is presented at the end of the script.

%% Closing and clearing the workspace

close all
clear all
clc

%% INPUT: Load data

% Input cell size & window size (in minutes)
cellsize=0.25;      % INPUT: Have to be the same with actual cell size
                    % Here it is 0.25 minutes (15 seconds)
cs=1/cellsize;    
winsize=8*cs;       % INPUT: 8 minutes = ~14.8 km
					% Window size can be modified but 
					% the number have to be even (6, 8, 10, etc.)

ss=cellsize*60;     % Cellsize in seconds
dd=cellsize/60;     % Cellsize in degrees

% INPUT: Gridded multibeam bathymetry with depth in METRES
% The resolution must comply with the determined 'cellsize'

data_input=load('Input_Bathymetry_15s.xyz'); 

lon=data_input(:,1);
lat=data_input(:,2);
depth=data_input(:,3);

%% Depth reshaped

l=length(find(lat==lat(1)));
lon_rs1=reshape(lon,l,[]);
lat_rs1=reshape(lat,l,[]);
depth_rs1=reshape(depth,l,[]);

%% Data borders for windowing

% Add top and bottom
lon_Ni=[lon_rs1(1,1)-dd*winsize:dd:lon_rs1(1,1)-dd]';
lon_Si=[lon_rs1(end,1)+dd:dd:lon_rs1(end,1)+dd*winsize]';

lat_NSi=lat_rs1(1,:);

for i=1:winsize
    for j=1:length(lon_rs1(1,:))
        depth_NS(i,j)=NaN;
        mask_NS(i,j)=NaN;
        lon_N(i,j)=lon_Ni(i);
        lon_S(i,j)=lon_Si(i);
        lat_NS(i,j)=lat_NSi(j);
    end
end
    
depth_input2=[depth_NS; depth_rs1; depth_NS];
lon_input2=[lon_N; lon_rs1; lon_S];
lat_input2=[lat_NS; lat_rs1; lat_NS];

% Add left and right
lat_Wi=fliplr([lat_input2(1,1)+dd:dd:lat_input2(1,1)+dd*winsize]);
lat_Ei=fliplr([lat_input2(1,end)-dd*winsize:dd:lat_input2(1,end)-dd]);

lon_EWi=lon_input2(:,1);

for i=1:winsize
    for j=1:length(depth_input2(:,1))
        depth_EW(i,j)=NaN;
        mask_EW(i,j)=NaN;
        lat_E(i,j)=lat_Ei(i);
        lat_W(i,j)=lat_Wi(i);
        lon_EW(i,j)=lon_EWi(j);        
    end
end

depth_rs=[depth_EW' depth_input2 depth_EW'];
lon_rs=[lon_EW' lon_input2 lon_EW'];
lat_rs=[lat_W' lat_input2 lat_E'];

%% LoG filter: Mask building

depth_im=mat2gray(depth_rs);

h = fspecial('log',winsize+1,0.1);      
I2 = imfilter(depth_im,h);

for i=1:length(I2(:,1))
    for j=1:length(I2(1,:))
        if I2(i,j) > 0
            log_mask(i,j)=NaN;
        else
            log_mask(i,j)=0;
        end
    end
end

figure()

x0=0; y0=0;
width=1000; height=800;
set(gcf,'position',[x0,y0,width,height])

pcolor(lon_rs,lat_rs,log_mask); shading flat
colorbar; colormap(gca,'gray')
set(gca,'FontSize',16)
title ('LoG Mask')
ytickangle(90); axis equal
xlim([min(min(lon_rs)) max(max(lon_rs))])
ylim([min(min(lat_rs)) max(max(lat_rs))])
caxis([0 1])

%% Aspect and slope (azimuth and plunge) computation

gridrv=[60*60/ss 0 0];
[aspects,slope,gradN,gradE]=gradientm(depth_rs,gridrv);

for i=1:length(aspects(:,1))
    for j=1:length(aspects(1,:))
        if aspects(i,j) >= 0 && aspects(i,j) <= 270
            aspect(i,j)=aspects(i,j)+90;
        else
            aspect(i,j)=aspects(i,j)-270;
        end
    end
end

%% Depth, aspect, and slope visualisation

figure()

x0=0; y0=0;
width=1000; height=800;
set(gcf,'position',[x0,y0,width,height])

pcolor(lon_rs,lat_rs,depth_rs/1000); shading flat
colorbar; colormap(gca,'jet')
set(gca,'FontSize',16)
title ('Depth (km)')
ytickangle(90); axis equal
xlim([min(min(lon_rs)) max(max(lon_rs))])
ylim([min(min(lat_rs)) max(max(lat_rs))])

figure()

x0=0; y0=0;
width=1000; height=800;
set(gcf,'position',[x0,y0,width,height])

pcolor(lon_rs,lat_rs,aspect); shading flat
colorbar; colormap(gca,'hsv')
caxis([0 360]); set(gca,'FontSize',16)
title ('Aspect (\circ)')
ytickangle(90); axis equal

figure()

x0=0; y0=0;
width=1000; height=800;
set(gcf,'position',[x0,y0,width,height])

pcolor(lon_rs,lat_rs,slope); shading flat
colorbar; colormap(gca,'bone')
caxis([0 30]); set(gca,'FontSize',16)
title ('Slope (\circ)')
ytickangle(90); axis equal

%% Making huge matrix: All dip direction (azimuth)

for k=1:winsize:(length(aspect(:,1))-winsize)*winsize
    for l=1:winsize:(length(aspect(1,:))-winsize)*winsize
        mat_aspect(k:k+winsize-1,l:l+winsize-1)=aspect((((k-1)/winsize)+1):(((k-1)/winsize)+1)+winsize-1,(((l-1)/winsize)+1):(((l-1)/winsize)+1)+winsize-1);
    end
end

%% Making huge matrix: All plunge (slope)

for k=1:winsize:(length(slope(:,1))-winsize)*winsize
    for l=1:winsize:(length(slope(1,:))-winsize)*winsize
        mat_slope(k:k+winsize-1,l:l+winsize-1)=slope((((k-1)/winsize)+1):(((k-1)/winsize)+1)+winsize-1,(((l-1)/winsize)+1):(((l-1)/winsize)+1)+winsize-1);
    end
end

%% Divide new matrix into windowed cells

nlin=floor(length(mat_aspect(:,1))/winsize);
ncol=floor(length(mat_aspect(1,:))/winsize);

vcell=ones(nlin,1)'*winsize;
hcell=ones(ncol,1)'*winsize;

cells=mat2cell(mat_aspect, vcell, hcell);

cells_slope=mat2cell(mat_slope, vcell, hcell);

%% Calculate eccentricity in each cell

for k=1:length(cells(:,1))
    for l=1:length(cells(1,:))
        [m,n]=size(cells{k,l});
        aspr=cells{k,l};
        slpr=cells_slope{k,l};
        
        % xyz components
        xc=sind(aspr).*cosd(-slpr);
        yc=cosd(aspr).*cosd(-slpr);
        zc=sind(-slpr);
        
        xr(k,l)=nansum(nansum(sind(aspr).*cosd(-slpr)));
        yr(k,l)=nansum(nansum(cosd(aspr).*cosd(-slpr)));
        zr(k,l)=nansum(nansum(sind(-slpr)));
        
        % R: resultant
        R(k,l)=sqrt(xr(k,l).^2+yr(k,l).^2+zr(k,l).^2);
        R_m(k,l)=R(k,l)/length(aspr(~isnan(aspr)));;
        
        % Lon, lat
        lon_win(k,l)=lon_rs(k,l)+(((winsize/cs)/2)/60);
        lat_win(k,l)=lat_rs(k,l)-(((winsize/cs)/2)/60);
        
        % Theta: plunge
        if xr(k,l) >= 0 && yr(k,l) >= 0
            theta_r(k,l)=atan(xr(k,l)/yr(k,l));
        elseif xr(k,l) < 0 && yr(k,l) >=0
            theta_r(k,l)=atan(xr(k,l)/yr(k,l))+degtorad(360);
        else
            theta_r(k,l)=atan(xr(k,l)/yr(k,l))+degtorad(180);
        end
        
        theta_d(k,l)=radtodeg(theta_r(k,l));
        
        % X-Count and Y-Count
        Xc(k,l)=R_m(k,l).*(sind(theta_d(k,l)));
        Yc(k,l)=R_m(k,l).*(cosd(theta_d(k,l)));
        
        % Eigen components
        num_x=length((cells{k,l}(~isnan(cells{k,l}))));
        
        if num_x > 0
            num_x = num_x;
        else
            num_x = 1;
        end
        
		% The B matrix
        Bs=[nansum(nansum(xc.^2)) nansum(nansum(xc.*yc)) nansum(nansum(xc.*zc));
            nansum(nansum(yc.*xc)) nansum(nansum(yc.^2)) nansum(nansum(yc.*zc));
            nansum(nansum(zc.*xc)) nansum(nansum(zc.*yc)) nansum(nansum(zc.^2))];

        B=Bs/num_x;

        D=eig(B);
        
        % Eccentricity
        
        a(k,l)=D(2);	% Semi-major axis
        b(k,l)=D(3);	% Semi-minor axis
        
        e(k,l)=sqrt(1-a(k,l).^2/b(k,l).^2);
    end
end

%% Make weight matrix (1-sin(slope)) and mask in the same size as eccentricity

for k=1:(length(depth_rs(:,1))-winsize)
    for l=1:(length(depth_rs(1,:))-winsize)
        sin_slope1(k,l)=sind((-slope(k+winsize/2,l+winsize/2)))/-(sind(max(max(slope))));
        weight(k,l)=1-sin_slope1(k,l);
        mask_filt(k,l)=log_mask(k+winsize/2,l+winsize/2);
    end
end

%% SWE: Slope weighted eccentricity

SWE=e.*weight;

%% Mask for SWE

for i=1:length(SWE(:,1))
    for j=1:length(SWE(1,:))
    if mask_filt(i,j) >= 0
        SWE_filt(i,j)=SWE(i,j);
    else
        SWE_filt(i,j)=NaN;
    end
    end
end

%% Plot e

figure()

x0=0; y0=0;
width=1000; height=800;
set(gcf,'position',[x0,y0,width,height])

pcolor(lon_win,lat_win,e); shading flat
colorbar
colormap(gca,'jet')
caxis([0.5 1])
xlim([min(min(lon_rs)) max(max(lon_rs))])
ylim([min(min(lat_rs)) max(max(lat_rs))])
set(gca,'FontSize',16); ytickangle(90)
title ('Eccentricity'); axis equal % of > 20 slope (\circ)')

%% Plot weight: 1-sin(slope)

figure()

x0=0; y0=0;
width=1000; height=800;
set(gcf,'position',[x0,y0,width,height])

pcolor(lon_win,lat_win,weight); shading flat
colorbar
colormap(gca,'jet')
caxis([0.5 1])
xlim([min(min(lon_rs)) max(max(lon_rs))])
ylim([min(min(lat_rs)) max(max(lat_rs))])
set(gca,'FontSize',16); ytickangle(90)
title ('Weight = 1-sin(slope)'); axis equal % of > 20 slope (\circ)')

%% Plot e * slope (SWE)

figure()

x0=0; y0=0;
width=1000; height=800;
set(gcf,'position',[x0,y0,width,height])

pcolor(lon_win,lat_win,SWE); shading flat
colorbar
colormap(gca,'jet')
caxis([0.5 1])
xlim([min(min(lon_rs)) max(max(lon_rs))])
ylim([min(min(lat_rs)) max(max(lat_rs))])
set(gca,'FontSize',16); ytickangle(90)
title ('Slope-weighted eccentricity (SWE)'); axis equal % of > 20 slope (\circ)')

%% Plot e * weight in masked area

figure()

x0=0; y0=0;
width=1000; height=800;
set(gcf,'position',[x0,y0,width,height])

cells_count=ncol*nlin;

pcolor(lon_win,lat_win,SWE_filt); shading flat
colorbar
colormap(gca,jet)
caxis([0.5 1])
xlim([min(min(lon_rs)) max(max(lon_rs))])
ylim([min(min(lat_rs)) max(max(lat_rs))])
set(gca,'FontSize',16); ytickangle(90)
title ('Masked SWE'); axis equal % of > 20 slope (\circ)')

%% WRITING DOWN

len=nlin*ncol;
e_rs=reshape(e,[1,len]);
weight_rs=reshape(weight,[1,len]);
SWE=e.*weight;
SWE_rs=reshape(SWE,[1,len]);
SWE_filt_rs=reshape(SWE_filt,[1,len]);

for i=1:nlin
    for j=1:ncol
        lon_win_sum(i,j)=lon_win(i,j);
        lat_win_sum(i,j)=lat_win(i,j);
    end
end

lon_win_sum_rs=reshape(lon_win_sum,[1,len]);
lat_win_sum_rs=reshape(lat_win_sum,[1,len]);

%% Eccentricity

fid1=fopen('Output_eccentricity.xyz', 'w');
for i=1:length(e_rs)
    fprintf(fid1, '%f %f %f\n',lon_win_sum_rs(i), lat_win_sum_rs(i), e_rs(i));
end

%% Weight: 1-sin(slope)

fid1=fopen('Output_weight.xyz', 'w');
for i=1:length(weight_rs)
    fprintf(fid1, '%f %f %f\n',lon_win_sum_rs(i), lat_win_sum_rs(i), weight_rs(i));
end

%% SWE: Slope-weighted eccentricity

fid1=fopen('Output_SWE.xyz', 'w');
for i=1:length(SWE_rs)
    fprintf(fid1, '%f %f %f\n',lon_win_sum_rs(i), lat_win_sum_rs(i), SWE_rs(i));
end

%% Masked SWE

fid1=fopen('Output_SWE_masked.xyz', 'w');
for i=1:length(SWE_filt_rs)
    fprintf(fid1, '%f %f %f\n',lon_win_sum_rs(i), lat_win_sum_rs(i), SWE_filt_rs(i));
end

%% XYZ to GRD

% The resulting *.xyz data can be converted into *.grd using the xyz2grd 
% function in GMT (http://gmt.soest.hawaii.edu/doc/5.3.2/xyz2grd.html)

% Use:
% gmtinfo *.xyz
% to obtain -Rxmin/xmax/ymin/ymax

% Then use:
% xyz2grd *.xyz -G*.grd -I15s -Rxmin/xmax/ymin/ymax
