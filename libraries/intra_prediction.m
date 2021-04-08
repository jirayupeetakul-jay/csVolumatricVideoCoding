function [SAD_candidate y_prediction, y_buffer_left, y_buffer_up, y_buffer_dc] = intra_prediction(y, phi, rows, columns, sub_pixels, m, n, y_buffer_left, y_buffer_up, y_buffer_cp, y_buffer_dc)
    %___SUM OF ABSOLUTE DIFFERENCE___
    SAD_y_left_mode   = sum(imabsdiff(y,y_buffer_left));
    SAD_y_up_mode     = sum(imabsdiff(y,y_buffer_up(:,columns)));
    SAD_y_cp_mode     = sum(imabsdiff(y,y_buffer_cp));
    SAD_y_dc_mode     = sum(imabsdiff(y,y_buffer_dc));
    SAD_candidate     = min([SAD_y_left_mode SAD_y_dc_mode SAD_y_cp_mode SAD_y_up_mode]);
    
    if(SAD_candidate == SAD_y_left_mode)
        y_prediction   = y_buffer_left;
    elseif(SAD_candidate == SAD_y_up_mode)
        y_prediction   = y_buffer_up(:,columns);
    elseif(SAD_candidate == SAD_y_dc_mode)
        y_prediction   = y_buffer_dc;
    else
        y_prediction = y_buffer_cp;
    end

    coefficient_matrix = sum(phi')';

    %___NEXT ROUND INTRA PREDICTION___
    switch n
        case 16
            if (rows <= 128) && (columns == 1)
                y_buffer_up(:,columns) = (((y(1)-y(2))/((sub_pixels*sub_pixels)/2)*coefficient_matrix));
                y_buffer_dc         = (y(1)/((sub_pixels*sub_pixels))*coefficient_matrix);
                y_buffer_left       = zeros(m, 1);
            elseif (rows == 1) && (columns <= 128)
                y_buffer_left       = ((y(1)-y(8))/((sub_pixels*sub_pixels)/2)*coefficient_matrix);
                y_buffer_dc         = (y(1)/((sub_pixels*sub_pixels))*coefficient_matrix);
                y_buffer_up(:,columns) = zeros(m, 1);
            else    
                y_buffer_left       = ((y(1)-y(8))/((sub_pixels*sub_pixels)/2)*coefficient_matrix);
                y_buffer_up(:,columns) = (((y(1)-y(2))/((sub_pixels*sub_pixels)/2)*coefficient_matrix));
                y_buffer_dc         = (y(1)/((sub_pixels*sub_pixels))*coefficient_matrix);
            end
            
        case 128
            if (rows <= 64) && (columns == 1)
                y_buffer_left       = zeros(m, 1) + (y(1)-y(2))/(2*sub_pixels);
                y_buffer_up(:,rows) = zeros(m, 1);
            elseif (rows == 1) && (columns <= 64) 
                y_buffer_up(:,rows) = zeros(m, 1) + (y(1)-y(m-3))/(2*sub_pixels);
                y_buffer_left       = zeros(m, 1);
            else
                y_buffer_left       = zeros(m, 1) + (y(1)-y(2))/(2*sub_pixels);
                y_buffer_up(:,rows) = zeros(m, 1) + (y(1)-y(m-3))/(2*sub_pixels);
            end 
            
        case 256
            %___NEXT ROUND INTRA PREDICTION___
            if (rows <= 2) && (columns == 1)
                y_buffer_left       = zeros(m, 1) + (y(1)-y(2))/(2*sub_pixels);
                y_buffer_up(:,rows) = zeros(m, 1);
            elseif (rows == 1) && (columns <= 256) 
                y_buffer_up(:,rows) = zeros(m, 1) + (y(1)-y(m-3))/(2*sub_pixels);
                y_buffer_left       = zeros(m, 1);
            else
                y_buffer_left       = zeros(m, 1) + (y(1)-y(2))/(2*sub_pixels);
                y_buffer_up(:,rows) = zeros(m, 1) + (y(1)-y(m-3))/(2*sub_pixels);
            end  
    end
end