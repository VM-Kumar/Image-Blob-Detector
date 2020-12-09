%TeamMembers
%MUTHU KUMAR VENKATESH:mvenkat3@ncsu.edu - 200318930
%SANJANA BANERJEE:sbaner24@ncsu.edu - 200316758
%KARTHIKA VADIVEL:kvadive@ncsu.edu -200322198
clc;
clear all;
close all;

%Reading the image for detection of blobs
image=imread('D:\DIS\PROJECT4\Project04_Final\Project04\TestImages\butterfly.jpg');
%to obtain program run time
tic;
%converting image to greyscale and double
image=rgb2gray(image);
image=im2double(image);

%initializing scale,sigma,scale space
[h w]= size(image);
scale_space = zeros(h,w,12);
k=1.32;
sigma=2;
       
   
%obtaining scale space:
%1.generating LoG kernel
%2.filtering the image using Log
%3.save square of laplacian
for i=1 : 12
           
filter=coeff( 2*ceil(sigma*3)+1, sigma);
filter=(sigma.^2)*filter;  
image_new=frequency_filter(image,filter);          

image_new=image_new.^2;  

if i==1
scale_space=image_new;
else
scale_space=cat(3,scale_space,image_new);
end
sigma=k*sigma;
           
           
end
   
   
%non maxima supression
%calculating maxima in each layer wrt neighbours
mx_new=[];
for i=1:12
img_layer = scale_space(:,:,i);        
[l m]=size(img_layer);
A=zeros(l+2,m+2);
A(2:l+1,2:m+1)=img_layer;
mx =zeros(size(img_layer));  
    for i= 2:l+1
        for j= 2:m+1
            mx(i-1,j-1)= max([A(i,j),A(i,j+1),A(i,j-1),A(i-1,j),A(i+1,j),A(i-1,j-1),A(i+1,j+1),A(i-1,j+1),A(i+1,j-1)]);
        end
    end
   
if i==1
mx_new=mx;
else
mx_new=cat(3,mx_new,mx);
end

end


%obtaining maximum among all the scale spaces
nms_3d=max(mx_new,[],3);
nms_3d= (nms_3d==mx_new).*mx_new;


%calculating the radius correspoiding to each sigma   
sigma=2;  
for i=1: 12
        radius(1,i)=1.414 * sigma;
               
        sigma=k*sigma;
end
   

%obtaining the centre coordinates, radius corresponding to the blobs
threshold=.007;  
for i=1: 12
img_layer = scale_space(:,:,i);
img_layer = (img_layer==nms_3d(:,:,i))&(img_layer>threshold);
[r,c] = find(img_layer);
    if i==1
    r_new=r;
    c_new=c;
    radius1=radius(i);        
    radius1=repmat(radius(i),size(r,1),1);
    else
    radius2=repmat(radius(i),size(r,1),1);
    radius1=cat(1,radius1,radius2);
    r_new=cat(1,r_new,r);
    c_new=cat(1,c_new,c);
    end
end

 
%display the elapsed time
toc;  
%display all circles corresponding to the blobs
display_circles(image, c_new, r_new, radius1);


%function to display the circle on the image
%I=input image
%Cx,Cy coordinates of the centre of circle
%radius=radius of the circle
function display_circles(I, cx, cy, radius1)

imshow(I);
hold on;
theta = 0:0.1:(2*pi+0.1);
cx1 = cx(:,ones(size(theta)));
cy1 = cy(:,ones(size(theta)));
rad = radius1(:,ones(size(theta)));
theta = theta(ones(size(cx1,1),1),:);
X = cx1+cos(theta).*rad;
Y = cy1+sin(theta).*rad;
line(X', Y', 'Color','y', 'LineWidth',1.5);
title(sprintf('%d circles', size(cx,1)));
end


%function to calculate filter kernal
%size=dimension of kernal
%s=sigma
function A = coeff(size,s)
L=(size/2)-0.5;
    for i=-L:L
        for j=-L:L
            A(i+L+1,j+L+1)=(-1/(pi*s^4))*(1-((i^2+j^2)/(2*s^2)))*exp(-(i^2+j^2)/(2*s^2));
        end
    end
end


%function to perform image filtering/2d convolution
%image=input image
%kernel=kernal to convolve with
 function outt=frequency_filter(image,kernel)
A1=im2double(image);
sz=size(A1);
A2=zeros(size(A1));
pad=[sz(1)-size(kernel,1), sz(2)- size(kernel,2)];
rowpad=floor((pad(1)+1)/2);
colpad=floor((pad(2)+1)/2);
A2(rowpad+1 : rowpad+size(kernel,1) , colpad+1 : colpad+size(kernel,2))=kernel;
kernel = ifftshift(A2);
x=D_2_fft(A1);
y=D_2_fft(kernel);
z=x.*y;
s=size(z);
outt=conj(D_2_fft(conj(z)))/(s(1)*s(2));
outt=outt*255;

%function to perform 2D dft using 1D dft
function img_fft2=D_2_fft(image)
L=size(image);
img_fft1=[];
    for i=1:L(1)
        row_ind=(image(i,:));
        partial_fft1=fft(row_ind);
        img_fft1=[img_fft1;partial_fft1];
    end
img_fft2=[];
    for j=1:L(2)
        col_ind=(img_fft1(:,j));
        partial_fft2=fft(col_ind);
        img_fft2=[img_fft2 partial_fft2];
    end
img_fft2;
fft2(image);
end
end