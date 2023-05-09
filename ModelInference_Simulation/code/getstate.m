function[ans]= getstate(x);
    if string(x)=='A' ans = 'O';
    elseif string(x)=='H' ans = 'S';
    elseif string(x)=='R' ans = 'S';
    elseif string(x)=='S' ans = 'S';
    elseif string(x)=='G' ans = 'A';
    else ans = 'A';
    end
end
