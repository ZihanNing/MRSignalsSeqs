function [M,Mxy] = epg2mag(Q)

M_mp = epg_FZ2spins(Q);
[M,Mxy] = epg_combM_convfrmQ(M_mp);