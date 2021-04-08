%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%MODE FUNCTIONS DEFINITIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef intra_prediction_9_modes
   methods (Static)
        function out=vertical_replication(T,N)
            for i=1:N
               for j=1:N
                   out(i,j)=T(i); % Vertical Replication
               end
            end
        end
        function out=horizonatal_replication(L,N)
            for i=1:N
                for j=1:N
                   out(i,j)=L(i);  % Horizonatal Replication
                end
            end
        end
        function out=mean_dc(L,T,N)
            for i=1:N
                for j=1:N
                    out(i,j)=round(mean([L(1:N) T(1:N) 4]));  % Mean / DC
                end
            end
        end
        function out=diagonal_down_left(T)
        %Diagonal Down Left
         a = (T(1) + 2*T(2) + T(3) + 2) / 4;
         b = (T(2) + 2*T(3) + T(4) + 2) / 4;
         c = (T(3) + 2*T(4) + T(5) + 2) / 4;
         d = (T(4) + 2*T(5) + T(6) + 2) / 4;
         e = (T(5) + 2*T(6) + T(7) + 2) / 4;
         f = (T(6) + 2*T(7) + T(8) + 2) / 4;
         g = (T(7) + 3*T(8)      + 2) / 4;

         out(1,1)=a;out(1,2)=b;out(1,3)=c;out(1,4)=d;
         out(2,1)=b;out(2,2)=c;out(2,3)=d;out(2,4)=e;
         out(3,1)=c;out(3,2)=d;out(3,3)=e;out(3,4)=f;
         out(4,1)=d;out(4,2)=e;out(4,3)=f;out(4,4)=g;

        end
        function out=diagonal_down_left_big(L,T,LT,N)
        %Plane for N=8 & 16
         H_bar = 1* (T(9) - T(7)) + 2* (T(10) - T(6))+ 3*(T(11) - T(5))+4*(T(12) - T(4))+ 5*(T(13)- T(3))+ 6*(T(14) - T(2))+ 7*(T(15)- T(1))+ 8*(T(16)- LT);

         V_bar = 1* (L(9) - L(6)) + 2* (L(10) - L(6))+ 3*(L(11) - L(5))+4*(L(12) - L(4))+ 5*(L(13)- L(3))+ 6*(L(14) - L(2))+ 7*(L(15)- L(1))+ 8*(L(16)- LT);

         H = (5*H_bar + 32) / 64;
         V = (5*V_bar + 32) / 64;

         a = 16 * (L(16) + T(16) + 1) - 7*(V+H);
         for j = 1: N
           for i = 1: N
             b = a + V * j + H * i;
             out(i,j) = (b/32);
           end
         end
        end
        function out=diagonal_down_right(L,T,LT)
        %diagonal Down right
         a = (L(4) + 2*L(3) + L(2) + 2) / 4;
         b = (L(3) + 2*L(2) + L(1) + 2) / 4;
         c = (L(2) + 2*L(1) + LT + 2) / 4;
         d = (L(1) + 2*LT + T(1) + 2) / 4;
         e = (LT + 2*T(1) + T(2) + 2) / 4;
         f = (T(1) + 2*T(2) + T(3) + 2) / 4;
         g = (T(1) + 2*T(3) + T(4) + 2) / 4;

         out(1,1)=d;out(1,2)=e;out(1,3)=f;out(1,4)=g;
         out(2,1)=c;out(2,2)=d;out(2,3)=e;out(2,4)=f;
         out(3,1)=b;out(3,2)=c;out(3,3)=d;out(3,4)=e;
         out(4,1)=a;out(4,2)=b;out(4,3)=c;out(4,4)=d;

        end

        function out=vertical_right(L,T,LT)
        %vertical right
         a = (LT + T(1) + 1) / 2;
         b = (T(1) + T(2) + 1) / 2;
         c = (T(2) + T(3) + 1) / 2;
         d = (T(3) + T(4) + 1) / 2;
         e = (L(1) + 2*LT + T(1) + 2) / 4;
         f = (LT + 2*T(1) + T(2) + 2) / 4;
         g = (T(1) + 2*T(2) + T(3) + 2) / 4;
         h = (T(2) + 2*T(3) + T(4) + 2) / 4;
         i = (LT + 2*L(1) + L(2) + 2) / 4;
         j = (L(1) + 2*L(2) + L(3) + 2) / 4;

         out(1,1)=a;out(1,2)=b;out(1,3)=c;out(1,4)=d;
         out(2,1)=e;out(2,2)=f;out(2,3)=g;out(2,4)=h;
         out(3,1)=i;out(3,2)=a;out(3,3)=b;out(3,4)=c;
         out(4,1)=j;out(4,2)=e;out(4,3)=f;out(4,4)=g;

        end

        function out=horizontal_down(L,T,LT)
        %Horizontal Down
         a = (LT + L(1) + 1) / 2;
         b = (L(1) + 2*LT + T(1) + 2) / 4;
         c = (LT + 2*T(1) + T(2) + 2) / 4;
         d = (T(1) + 2*T(2) + T(3) + 2) / 4;
         e = (L(1) + L(2) + 1) / 2;
         f = (LT + 2*L(1) + L(2) + 2) / 4;
         g = (L(2) + L(3) + 1) / 2;
         h = (L(1) + 2*L(2) + L(3) + 2) / 4;
         i = (L(3) + L(4) + 1) / 2;
         j = (L(2) + 2*L(3) + L(4) + 2) / 4;

         out(1,1)=a;out(1,2)=b;out(1,3)=c;out(1,4)=d;
         out(2,1)=e;out(2,2)=f;out(2,3)=a;out(2,4)=b;
         out(3,1)=g;out(3,2)=h;out(3,3)=e;out(3,4)=f;
         out(4,1)=i;out(4,2)=j;out(4,3)=g;out(4,4)=h;

        end
        function out=vertical_left(T)
        %Vertical Left
         a = (T(1) + T(2) + 1) / 2;
         b = (T(2) + T(3) + 1) / 2;
         c = (T(3) + T(4) + 1) / 2;
         d = (T(4) + T(5) + 1) / 2;
         e = (T(5) + T(6) + 1) / 2;
         f = (T(1) + 2*T(2) + T(3) + 2) / 4;
         g = (T(2) + 2*T(3) + T(4) + 2) / 4;
         h = (T(3) + 2*T(4) + T(5) + 2) / 4;
         i = (T(4) + 2*T(5) + T(6) + 2) / 4;
         j = (T(5) + 2*T(6) + T(7) + 2) / 4;

         out(1,1)=a;out(1,2)=b;out(1,3)=c;out(1,4)=d;
         out(2,1)=f;out(2,2)=g;out(2,3)=h;out(2,4)=i;
         out(3,1)=b;out(3,2)=c;out(3,3)=d;out(3,4)=e;
         out(4,1)=g;out(4,2)=h;out(4,3)=i;out(4,4)=j;

        end
        function out=horizontal_up(L)
        %Horizontal UP
         a = (L(1) + L(2) + 1) / 2;
         b = (L(1) + 2*L(2) + L(3) + 2) / 4;
         c = (L(2) + L(3) + 1) / 2;
         d = (L(2) + 2*L(3) + L(4) + 2) / 4;
         e = (L(3) + L(4) + 1) / 2;
         f = (L(3) + 3*L(4)      + 2) / 4;
         g = L(4);

         out(1,1)=a;out(1,2)=b;out(1,3)=c;out(1,4)=d;
         out(2,1)=c;out(2,2)=d;out(2,3)=e;out(2,4)=f;
         out(3,1)=e;out(3,2)=f;out(3,3)=g;out(3,4)=g;
         out(4,1)=g;out(4,2)=g;out(4,3)=g;out(4,4)=g;

        end
        
        function intra_prediction_candidate=sum_of_absolute_difference(c1,c2,c3,c4,c5,c6,c7,c8,c9,sub_frame_temp)
            SAD(1)=sum(sum(abs(c1-sub_frame_temp)));
            SAD(2)=sum(sum(abs(c2-sub_frame_temp)));
            SAD(3)=sum(sum(abs(c3-sub_frame_temp)));
            SAD(4)=sum(sum(abs(c4-sub_frame_temp)));
            SAD(5)=sum(sum(abs(c5-sub_frame_temp)));
            SAD(6)=sum(sum(abs(c6-sub_frame_temp)));
            SAD(7)=sum(sum(abs(c7-sub_frame_temp)));
            SAD(8)=sum(sum(abs(c8-sub_frame_temp)));
            SAD(9)=sum(sum(abs(c9-sub_frame_temp)));
            
            % Selection based on Min SAD
            [~,min_SAD]=min(SAD);
            switch min_SAD
                case 1
                   fprintf('Vertical Replication \n');
                    intra_prediction_candidate=c1;
                    intra_mode='Mode 1';
                case 2
                   fprintf('Horizonatal Replication \n');
                    intra_prediction_candidate=c2;
                    intra_mode='Mode 2';
                case 3
                    fprintf('Mean/DC \n');
                    intra_prediction_candidate=c3;
                    intra_mode='Mode 3';
                case 4
                    fprintf('Diagonal Down-Left \n');
                    intra_prediction_candidate=c4;
                    intra_mode='Mode 4';
                case 5
                    fprintf('Diagonal Down-Right \n');
                    intra_prediction_candidate=c5;
                    intra_mode='Mode 5';
                case 6
                    fprintf('Vertical Right \n');
                    intra_prediction_candidate=c6;
                    intra_mode='Mode 6';
                case 7
                    fprintf('Horizontal Down \n');
                    intra_prediction_candidate=c7;
                    intra_mode='Mode 7';
                case 8
                    fprintf('Vertical Left \n');
                    intra_prediction_candidate=c8;
                    intra_mode='Mode 8';
                case 9
                    fprintf('Horizontal Up \n');
                    intra_prediction_candidate=c9;
                    intra_mode='Mode 9';
            end
        end
   end
end