function [M_com,Mxy] = epg_combM_convfrmQ(M)
    Mx = sum(M(1,:));
    My = sum(M(2,:));
    Mz = sum(M(3,:));
    M_com = [Mx,My,Mz]';

    Mxy = Mx + 1i * My;
%     signal_magnitude = abs(Mxy);   
%     signal_phase = angle(Mxy);     

end