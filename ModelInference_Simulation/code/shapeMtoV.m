function[Avec] = shapeMtoV(A,ii,jj,kk)
    Avec = zeros(1,ii*jj*kk);
    for id=1:ii*jj*kk
        id_tmp = id;
        i1 = floor((id_tmp-1)/jj/kk);
        id_tmp = id_tmp - jj*kk*i1;
        i2 = floor((id_tmp-1)/kk);
        id_tmp = id_tmp - kk*i2;
        i3 = id_tmp;
        i1 = i1+1;
        i2 = i2+1;
        Avec(id) = A(i1,i2,i3);
    end
end


