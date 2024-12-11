function [x,y,rev,cTickets,cSeats2,cSeats3] = calc_decision(fp1,fp2j,fp3j,r1,r2j,r3j,dec,x,y,t,n,cTickets,cSeats2,cSeats3)
    % Let the customer decide
    if dec(t,n) > (fp1 + sum(fp2j) + sum(fp3j))
        % Leave without a ticket
        cTickets(1) = cTickets(1)+1;
        rev = 0;
        return
    else
        concat = [fp1; fp2j(:); fp3j(:)];
        cumsumConcat = cumsum(concat);
        lower = dec(t,n) < cumsumConcat;
        ind = find(lower,1,'first');
        if ind == 1
            % Buy product 1 (without reservation)
            y = y - 1;
            rev = r1;
            cTickets(2) = cTickets(2)+1;
        elseif ind <= numel(fp2j) + 1
            % Buy product 2 at seat ind
            x(ind - 1) = 0;
            rev = r2j(ind-1);
            y = y - 1;
            cTickets(3) = cTickets(3)+1;
            cSeats2(ind-1) = cSeats2(ind-1)+1; 
        else
            % Buy product 3 at seat ind and its adjacent one
            ind = ind - numel(fp2j) - 1;
            x(ind) = 0;
            x(ind - (-1)^ind) = 0;
            rev = r3j(ind);
            y = y - 2;
            cSeats3(ind) = cSeats3(ind)+1; 
            cSeats3(ind - (-1)^ind) = cSeats3(ind - (-1)^ind)+1;
            cTickets(4) = cTickets(4)+1;
        end
    end
end
