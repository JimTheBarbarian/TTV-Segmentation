function [csf_dice, gm_dice, wm_dice] = compute_brain_dice(segment_u, gt)
   
    %get ground truth segmentation for each brain portion
    csf_m1 = double(gt== 48);
    white_matter_m1 = double(gt==154);
    grey_matter_m1 = double(gt==106); 

    %map segmentation result to a clustering result
    [~,threshold_u] = max(segment_u,[],3);
    
    %compute each dice
    csf_dice = max([dice(double(threshold_u==1), csf_m1), dice(double(threshold_u==2), csf_m1), dice(double(threshold_u==3), csf_m1), dice(double(threshold_u==4), csf_m1)]);
    wm_dice = max([dice(double(threshold_u==1), white_matter_m1), dice(double(threshold_u==2), white_matter_m1), dice(double(threshold_u==3), white_matter_m1), dice(double(threshold_u==4), white_matter_m1)]);
    gm_dice = max([dice(double(threshold_u==1), grey_matter_m1), dice(double(threshold_u==2), grey_matter_m1), dice(double(threshold_u==3), grey_matter_m1), dice(double(threshold_u==4), grey_matter_m1)]);

    

end