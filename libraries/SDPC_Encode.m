

function [i y_q y_index]= SDPC_Encode(y,q,DC_Measure)


i = zeros(size(y));
y_q = zeros(size(y));

y_hat = DC_Measure;

y_index = zeros(size(y,2),1);
L_norm = 1; % 1 or 2
count_l = 0; % left   index = 0
count_u = 0; % up     index = -1
count_ul = 0;% upleft index = 2
count_d = 0; % dc     index = 1
coef = zeros(size(y,2),1);

for j = 1:size(y,2)
    if j == 1
        
        d = y(:,j) - y_hat;
        coef(j) = dot(y(:,j),y_hat)/(norm(y(:,j),2)*norm(y_hat,2)+eps)
        i(:,j) = round(d/q);
        d_hat = i(:,j)*q;
        y_q(:,j) = y_hat + d_hat;
        y_index(j) = 0;
        
    elseif ceil(j/32) == 1
        
        value_l = norm(y(:,j)-DC_Measure,L_norm);
        value_u = norm(y(:,j)-y_q(:,j-1),L_norm);
        value_d = norm(y(:,j)-((DC_Measure+y_q(:,j-1))/2),L_norm);
        [min_value min_index] = min([value_u,value_l,value_d]);
        
        switch min_index
            case 1
                y_hat = y_q(:,j-1);
                count_u = count_u + 1;
                y_index(j) = -1;
            case 2
                y_hat = DC_Measure;
                count_l = count_l + 1;
                y_index(j) = 0;
            case 3
                y_hat = (DC_Measure + y_q(:,j-1))/2;
                count_d = count_d + 1;
                y_index(j) = 1;
        end
        
        coef(j) = dot(y(:,j),y_hat)/(norm(y(:,j),2)*norm(y_hat,2)+eps);
        d = y(:,j) - y_hat;
        i(:,j) = round(d/q);
        d_hat = i(:,j)*q;
        y_q(:,j) = y_hat + d_hat;
        
    elseif mod(j,32) == 1
        
        value_l = norm(y(:,j)-y_q(:,j-32),L_norm);
        value_u = norm(y(:,j)-DC_Measure,L_norm);
        value_d = norm(y(:,j)-((DC_Measure+y_q(:,j-32))/2),L_norm);
        [min_value min_index] = min([value_u,value_l,value_d]);
        
        switch min_index
            case 1
                y_hat = DC_Measure;
                count_u = count_u + 1;
                y_index(j) = -1;
            case 2
                y_hat = y_q(:,j-32);
                count_l = count_l + 1;
                y_index(j) = 0;
            case 3
                y_hat = (DC_Measure + y_q(:,j-32))/2;
                count_d = count_d + 1;
                y_index(j) = 1;
        end
        
        coef(j) = dot(y(:,j),y_hat)/(norm(y(:,j),2)*norm(y_hat,2)+eps);
        d = y(:,j) - y_hat;
        i(:,j) = round(d/q);
        d_hat = i(:,j)*q;
        y_q(:,j) = y_hat + d_hat;
    else
        
        %if norm(y(:,j)-y_q(:,j-1),2) <= norm(y(:,j)-y_q(:,j-32),2)
        value_l = norm(y(:,j)-y_q(:,j-32),L_norm);
        value_u = norm(y(:,j)-y_q(:,j-1),L_norm);
        value_ul = norm(y(:,j)-y_q(:,j-33),L_norm);
        value_d = norm(y(:,j)-((y_q(:,j-32)+y_q(:,j-1))/2),L_norm);
        [min_value min_index] = min([value_u,value_l,value_ul,value_d]);
        
        switch min_index
            case 1
                y_hat = y_q(:,j-1);
                count_u = count_u + 1;
                y_index(j) = -1;
            case 2
                y_hat = y_q(:,j-32);
                count_l = count_l + 1;
                y_index(j) = 0;
            case 3
                y_hat = y_q(:,j-33);
                count_ul = count_ul + 1;
                y_index(j) = 2;
            case 4
                y_hat = (y_q(:,j-32) + y_q(:,j-1))/2;
                count_d = count_d + 1;
                y_index(j) = 1;
                
        end
        
        coef(j) = dot(y(:,j),y_hat)/(norm(y(:,j),2)*norm(y_hat,2)+eps);
        d = y(:,j) - y_hat;
        i(:,j) = round(d/q);
        d_hat = i(:,j)*q;
        y_q(:,j) = y_hat + d_hat;
        
    end
    
end


