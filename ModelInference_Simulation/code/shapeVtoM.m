function[A] = shapeVtoM(Avec,ii,jj,kk)
    A = zeros([ii,jj,kk]);
    for id1 = 1:ii
        for id2=1:jj
            for id3=1:kk
                A(id1,id2,id3) = Avec((id1-1)*jj*kk + (id2-1)*jj + id3);
            end
        end
    end
end
