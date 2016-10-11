
[h,d] = read('Thalamus-maxprob-thr25-2mm.nii');

right = d;
left = d;

for i = 1:45
    
    left(i,:,:) = 0;
    
end

for i = 45:size(d,1)
    
    
    right(i,:,:) = 0;
    
end

write(h,left,'left.nii');
write(h,right,'right.nii');