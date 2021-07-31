function lv = hsslevel(tr)

switch length(tr)
    case 0
        error('error');
    case 1
        lv = 0;
    case 2
        lv = [1 0];
    case 3
        lv = [1 1 0];
    otherwise
        lv = 0;
        ch = child(tr); %sort(find(tr == length(tr)));
        n = length(tr);
        lc = ch{n}(1);
        while ~isempty(ch{lc})
            lc = ch{lc}(1);
        end
        if length(ch{n}) == 2
            lv1 = hsslevel(tr(lc:ch{n}(1)));
            lv2 = hsslevel(tr(ch{n}(1)+1:ch{n}(2))-ch{n}(1));
            lv = [lv1+1 lv2+1 0];
        elseif length(ch{n}) == 1
            lv = [hsslevel(tr(lc:ch{n}(1)))+1 0];
        end
end