function until = StraightUntil(X,Y,start_idx,dir)

par1 = numel(X)/20;
reg = 3;

avgs = zeros(1,round(numel(X)/reg));

if dir == 1
    interest = Y(start_idx:end);
    accp_data = 0;
    i = 1;
    stop = 0;

    while i <= numel(interest)-3 && ~isnan(accp_data) && stop == 0
        if mod(i,reg) == 0
            avgs(1,i/reg) = (interest(i)+interest(i-1)+interest(i-2))/reg;
            if i > reg && avgs(i/reg) > avgs(i/reg - 1) && accp_data>0
                accp_data = accp_data - reg;
                stop = 1;
            end
        end
        if interest(i) < interest(i+1) && interest(i+1) < interest(i+2) && interest(i+2) < interest(i+3)  && stop ~= 1
            accp_data = accp_data - 2;
            i = numel(X);
        else
            accp_data = accp_data + 1;
            i = i + 1;
        end
    end

else
    interest = Y(1:start_idx);

    accp_data = 0;
    i = 1;
    stop = 0;

    while i <= numel(interest)-3 && ~isnan(accp_data) && stop == 0
        if mod(i,reg) == 0
            avgs(1,i/reg) = (interest(end-(i-3))+interest(end-(i-2))+interest(end-(i-1))/reg);
            if i > reg && avgs(i/reg) < avgs(i/reg - 1) && accp_data>0
                accp_data = accp_data - reg;
                stop = 1;
            end
        end
        if interest(end-(i-1)) > interest(end-i) && interest(end-i) > interest(end-(i+1)) && interest(end-(i+1)) > interest(end-(i+2)) && stop ~= 1
            accp_data = accp_data -2;
            i = numel(X);
        else
            accp_data = accp_data + 1;
            i = i + 1;
        end
    end
end

if par1 <= accp_data && dir == 1
    until = start_idx + accp_data;
elseif par1 <= accp_data && dir == -1
    until = start_idx - accp_data;
else
    until = NaN;
end


end
