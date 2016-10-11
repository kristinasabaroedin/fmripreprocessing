function [mm_coords]=get_mm_coords(filename)
% This script will spit out all the coordinates of a mask in MNI mm coordinates. 
% 
% The only thing to watch for is whether you want neurological versus radiological orientation. 
% Just delete the line with x = -x to toggle between the two orientations. 
% Left brain is negative if you keep x = -x. Left brain is positive if you delete x = -x. 

% filename = ('R.nii'); % filename of image
voxdim = [2 2 2]; % voxel dimensions of input image (1 x 3 vector)
origin=[46,64,37]; % voxel coords of origin of MNI space i.e. x=0mm is voxel coordinate 46
orient = 1; % 1 for left being negative; 0 for left being positive.


[hdr,data]=read(filename); %read in mask filename
[x,y,z]=ind2sub(size(data),find(data)); %voxel coordinates
x=voxdim(1)*(x-origin(1)); %multiply by voxdim and subtract origin
y=voxdim(2)*(y-origin(2));
z=voxdim(3)*(z-origin(3));

if orient ==1
    x=-x; %neurological versus radiological
end

mm_coords = [x,y,z];
%x, y and z are now MNI coordinated expressed in mm 