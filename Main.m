%% 4F imaging system using fourier optics

clear all;
close all;
clc;

%% set the main parameters of all matrices and filters used

dim=512;      % the dimension of the array image or filter
center= dim/2 +1 ;
s= 3;
f=1/(dim*s); % frequency spacing 

%% reading the input image

img = imread('InputImage.jpg');
img = imresize(img,[dim,dim]);

%bin_object = rgb2gray(object);

%% creatig the filters

% 1.Horizontal Single Slit Aperture
H_SingleSlit = zeros(dim,dim);
v_single_slit_width   = 100;   % width pixels
H_SingleSlit_hight  = 20;    %  hight pixels

H_SingleSlit((center-H_SingleSlit_hight):(center+H_SingleSlit_hight),(center-v_single_slit_width):(center+v_single_slit_width)) =ones(2*H_SingleSlit_hight+1,2*v_single_slit_width+1);

% horizontal double slit 
horizontal_double_slit = zeros(dim,dim);
horizontal_double_slit_width   = 100;   
horizontal_double_slit_hight  = 20;    
horizontal_double_slit_gap     = 50;   
horizontal_double_slit(((center-horizontal_double_slit_gap/2)-horizontal_double_slit_hight):((center-horizontal_double_slit_gap/2)+horizontal_double_slit_hight),...
    (center-horizontal_double_slit_width):(center+horizontal_double_slit_width)) = ...
                                          ones(2*horizontal_double_slit_hight+1,2*horizontal_double_slit_width+1);
                                      
horizontal_double_slit(((center+horizontal_double_slit_gap/2)-horizontal_double_slit_hight):((center+horizontal_double_slit_gap/2)+horizontal_double_slit_hight),...
    (center-horizontal_double_slit_width):(center+horizontal_double_slit_width)) = ...
                                          ones(2*horizontal_double_slit_hight+1,2*horizontal_double_slit_width+1);

              

%% Vertical Single Slit Aperture
vertical_single_slit = zeros(dim,dim);
v_single_slit_width   = 20;   % pixels
vertical_single_slit_hight  = 100;    % pixels
vertical_single_slit((center-vertical_single_slit_hight):(center+vertical_single_slit_hight),(center-v_single_slit_width):(center+v_single_slit_width)) = ...
                                          ones(2*vertical_single_slit_hight+1,2*v_single_slit_width+1);
 
                                      
                                      
% Vertical Double Slit Aperture
vertical_double_slit = zeros(dim,dim);
v_single_slit_width   = 20;   % pixels
vertical_single_slit_hight  = 100;    % pixels
v_spacing     = 60;   % pixels
vertical_double_slit((center-vertical_single_slit_hight):(center+vertical_single_slit_hight),...
    ((center-v_spacing/2)-v_single_slit_width):((center-v_spacing/2)+v_single_slit_width)) = ...
                                          ones(2*vertical_single_slit_hight+1,2*v_single_slit_width+1);
                                      
vertical_double_slit((center-vertical_single_slit_hight):(center+vertical_single_slit_hight),...
    ((center+v_spacing/2)-v_single_slit_width):((center+v_spacing/2)+v_single_slit_width)) = ...
                                          ones(2*vertical_single_slit_hight+1,2*v_single_slit_width+1);

                                      
                                      
%% pinhole                                       
                                      
Xfreq = ((-dim/2):(dim/2-1))*f;
Yfreq = -Xfreq;
[FX,FY] = meshgrid(Xfreq,Yfreq);  %2-D arrays hold fx location and fy location of all points
freq_rad = sqrt(FX.^2 + FY.^2);
maxfreq = (dim/2-1)*f;

cutoff_freq1 = 0.1*maxfreq;
pinhole_filter = double(freq_rad <= cutoff_freq1); % pinhole
                                      
                                      
                                      
%% Fourier transform of the Image
FT_img= fftshift(fft2(fftshift(img(:,:,3))));
ftimg_H_SingleSlit = FT_img.*H_SingleSlit;
ftimg_horizontal_double_slit = FT_img.*horizontal_double_slit;

ftimg_vertical_single_slit = FT_img.*vertical_single_slit;
ftimg_vertical_double_slit = FT_img.*vertical_double_slit;

ft_pinhole = FT_img.*pinhole_filter;



% calculate inverse FT and the absolute value to plot it  
img1_H_SingleSlit = abs(fftshift(ifft2(fftshift(ftimg_H_SingleSlit))));
img2_horizontal_double_slit = abs(fftshift(ifft2(fftshift(ftimg_horizontal_double_slit))));

img3_vertical_single_slit = abs(fftshift(ifft2(fftshift(ftimg_vertical_single_slit))));
img4_vertical_double_slit = abs(fftshift(ifft2(fftshift(ftimg_vertical_double_slit))));
img_pinhole = abs(fftshift(ifft2(fftshift(ft_pinhole))));



%%    >>>>>> PLOTTING

% 1. horezontal single slit 

figure('Name', 'Horizontal Single slit');
%set(gcf, 'Units','Normalized','OuterPosition',[0 0 1 1]);
colormap('summer');

subplot(1,4,1);
imagesc(H_SingleSlit);axis('image');
title('Horozental single slit');

subplot(1,4,2);
imagesc(img1_H_SingleSlit);axis('image');
xlabel('fx, cycles/m');ylabel('fy, cycles/m');
colorbar('EastOutside');
title('horzental single sit');


% horizontal double slit 
figure('Name', 'Horizontal double slit');
colormap('summer');
subplot(1,4,1);
imagesc(horizontal_double_slit);axis('image');
title('Horizontal double slit');

subplot(1,4,2);
imagesc(img2_horizontal_double_slit);axis('image');
xlabel('fx, cycles/m');ylabel('fy, cycles/m');
colorbar('EastOutside');
title('horizontal double slit image');




% pinhole
figure('Name', 'pinhole');
colormap('summer');
subplot(1,4,1);
imagesc(pinhole_filter);axis('image');
title('pinhole filter');

subplot(1,4,2);
imagesc(img_pinhole);axis('image');
xlabel('fx, cycles/m');ylabel('fy, cycles/m');
colorbar('EastOutside');
title('pinhole Image');

                    
