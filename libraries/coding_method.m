function [y_predicted, mode] = coding_method(y, phi, i, j, sub_pixels, m, n, y_left, y_up, y_up_left, method)
    switch method
        case 'intra_prediction'
            %___SUM OF ABSOLUTE DIFFERENCE___
            SAD_y_left_mode   = norm(y_left - y);
            SAD_y_up_mode     = norm(y_up - y);
            SAD_y_cp_mode     = norm(y);
            SAD_y_dc_mode     = norm(round((y_left+y_up)/2) - y);
            [value, mode]     = min([SAD_y_cp_mode SAD_y_left_mode SAD_y_up_mode SAD_y_dc_mode]);
            
            if(mode == 1)
                y_predicted   = zeros(m,1);
            elseif(mode == 2)
                y_predicted   = y_left;
            elseif(mode == 3)
                y_predicted   = y_up;
            elseif(mode == 4)
                y_predicted   = round((y_up+y_left)/2);
                mode          = 9;
            else
                y_predicted   = zeros(m,1);
            end

            coefficient_matrix = sum(phi')'; % sum of transpose matrix

            %___NEXT ROUND INTRA PREDICTION___
            switch n
                case 256
                     y_buffer_left = floor((y(1)-y(3))/((sub_pixels*sub_pixels)/2))*coefficient_matrix;
                     y_buffer_up   = floor((y(1)-y(2))/((sub_pixels*sub_pixels)/2))*coefficient_matrix;
                     y_buffer_dc   = (y_buffer_left + y_buffer_up)/2;
            end
        case 'intra_frame_coding'
            SAD_y_cp_mode                = norm(y);
            SAD_y_left_mode              = norm(y_left - y); %A
            SAD_y_up_mode                = norm(y_up - y); %B
%             SAD_y_up_left_mode           = norm(y_up_left - y); %C
%             SAD_y_left_up_up_left_mode   = norm(round(y_left+y_up-y_up_left) - y); %A+B-C
%             SAD_y_left_up_up_left_2_mode = norm(round(y_left+round(y_up-y_up_left)/2) - y); %A+(B-C)/2
%             SAD_y_left_up_2_up_left_mode = norm(round((y_left+y_up)/2)-y_up_left - y); %(A+B)/2-C
%             SAD_y_up_left_up_left_2_mode = norm(round(y_up+((y_left-y_up_left)/2)) - y); %B+(A–C)/2 
            SAD_y_dc_mode                = norm(round((y_left+y_up)/2) - y);
%             [value, mode]                = min([SAD_y_cp_mode SAD_y_left_mode SAD_y_up_mode SAD_y_up_left_mode SAD_y_left_up_up_left_mode ...
%                                                 SAD_y_left_up_up_left_2_mode SAD_y_left_up_2_up_left_mode SAD_y_up_left_up_left_2_mode SAD_y_dc_mode]);
            [value, mode]                = min([SAD_y_cp_mode SAD_y_left_mode SAD_y_up_mode SAD_y_dc_mode]);
            if(mode == 1)
                y_predicted   = zeros(m,1);
            elseif(mode == 2)
                y_predicted   = y_left;
            elseif(mode == 3)
                y_predicted   = y_up;
%             elseif(mode == 4)
%                 y_predicted   = y_up_left;
%             elseif(mode == 5)
%                 y_predicted   = round(y_left+y_up-y_up_left);
%             elseif(mode == 6)
%                 y_predicted   = round(y_left+(y_up-y_up_left)/2);
%             elseif(mode == 7)
%                 y_predicted   = round((y_left+y_up)/2)-y_up_left;
%             elseif(mode == 8)
%                 y_predicted   = round(y_up+((y_left-y_up_left)/2));
            elseif(mode == 4)
                y_predicted   = round((y_up+y_left)/2);
            else
                y_predicted   = zeros(m,1);
            end
    end
end