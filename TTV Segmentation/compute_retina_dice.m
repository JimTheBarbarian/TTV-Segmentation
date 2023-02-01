function dice_val =  compute_retina_dice(u_segment, gt)

dice_val = dice(double(u_segment(:,:,1)>0.5), double(gt==255));
end