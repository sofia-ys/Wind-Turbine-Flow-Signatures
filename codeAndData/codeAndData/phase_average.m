load("V_velocity_matrix.mat")
AVG = zeros(146204, 24);
for i = 1:146204
    v_avg = 0;
    for ii = 1:24
        for iii = 0:99
            v_avg = v_avg + V_matrix.velocity(i, ii+iii*24);
        end
        v_avg = v_avg/100;
        AVG(i, ii)=v_avg;
    end
end
save('AVG.mat', 'AVG');
