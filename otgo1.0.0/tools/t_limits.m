function [t_lower, t_upper] = t_limits(ln_p1,lnc,min_max)
            % TLIMITS lower and upper limits denote the limits of the variable of line
            %  equtions t that intersects with the cell.
            % 
            % [TLOWER, TUPPER] = TLIMITS(CELL,LN) calculates lower and
            %  upper limit of t values of the line LN. This interval of t
            %  values denotes the portion of line LN that possibly intersets
            %  with the cell CELL.
            %
            % See also Cell.

            if lnc(1) > 0
                tx_lower = (min_max(1,1) - ln_p1(1)) / lnc(1);
                if tx_lower < 0 || round(lnc(1),6) == 0
                    tx_lower = 0;
                end
                tx_upper = (min_max(2,1) - ln_p1(1)) / lnc(1);
                if tx_upper < 0
                    tx_upper = Inf;
                end
            else
                tx_lower = (min_max(2,1) - ln_p1(1)) / lnc(1);
                if tx_lower < 0 || round(lnc(1),6) == 0
                    tx_lower = 0;
                end
                tx_upper = (min_max(1,1) - ln_p1(1)) / lnc(1);
                if tx_upper < 0
                    tx_upper = Inf;
                end
            end

            if lnc(2) > 0
                ty_lower = (min_max(1,2) - ln_p1(2)) / lnc(2);
                if ty_lower < 0 || round(lnc(2),6) == 0
                    ty_lower = 0;
                end
                ty_upper = (min_max(2,2) - ln_p1(2)) / lnc(2);
                if ty_upper < 0
                    ty_upper = Inf;
                end
            else
                ty_lower = (min_max(2,2) - ln_p1(2)) / lnc(2);
                if ty_lower < 0 || round(lnc(2),6) == 0
                    ty_lower = 0;
                end
                ty_upper = (min_max(1,2) - ln_p1(2)) / lnc(2);
                if ty_upper < 0
                    ty_upper = Inf;
                end
            end

            if lnc(3) > 0
                tz_lower = (min_max(1,3) - ln_p1(3)) / lnc(3);
                if tz_lower < 0 || round(lnc(3),6) == 0
                    tz_lower = 0;
                end
                tz_upper = (min_max(2,3) - ln_p1(3)) / lnc(3);
                if tz_upper < 0
                    tz_upper = Inf;
                end
            else
                tz_lower = (min_max(2,3) - ln_p1(3)) / lnc(3);
                if tz_lower < 0 || round(lnc(3),6) == 0
                    tz_lower = 0;
                end
                tz_upper = (min_max(1,3) - ln_p1(3)) / lnc(3);
                if tz_upper < 0
                    tz_upper = Inf;
                end
            end
            t_lower = max([tx_lower,ty_lower,tz_lower]);
            t_upper = min([tx_upper,ty_upper,tz_upper]);
        end