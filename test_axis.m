%% helper scrip demonstrating how the axis convention works
% This way x,y are alligned with the image and y grows downward (as
% expected). Thus z grows forward. The origin is on the image plane
% Author Fotios Logothetis fl302@cam.ac.uk
clear all
close all

nrows=108; %standard image size /10
ncols=192;
f=100;  %irrelevant
%lower case letter are image plane coordinates, uppercase are 3D word
[x,y] = meshgrid(1:ncols,1:nrows);

x=x-ncols/2;
y=y-nrows/2;

z=zeros(nrows,ncols);

Z = repmat(1:ncols,nrows,1); %get a non-constant z
X = x.*(f+Z)/f;
Y = y.*(f+Z)/f;

mask=ones(nrows,ncols);
mask(:,1:10)=0; %make a nice non symmetric triangle
for i=1:nrows
   mask(i,2*i:end)=0;    
end

figure;
imshow(mask);

Z(mask==0)=NaN;
z(mask==0)=NaN;

figure;
%show 3D shape
surf(X,Y,Z);

hold on; 
%show axis and origin
plot3(0,0,-f,'co');
for ii=1:5
	plot3( ii^2, 0,0,'r*');
    plot3( 0, ii^2, 0, 'g*');
	plot3( 0,0, ii^2, 'b*');
end
%show image plane
surf(x,y,z);

hold off;