perl /mnt/qjw/myperl/fst_pi_common.pl E_W_fst.windowed.weir.fst Eastern_high.windowpi.out.pi Western_low.windowpi.out.pi 10 |sort -k 1,1 -k 2,2n > E_W_fst_pi
perl /media/primates/wrf01/qjwperl/draw_Zfst_Pi.pl ./E_W_fst_pi E-W 0.05
