#0.5 Hz
#Amylin
#SERCA=100%
#python daisychain.py -dt 0.1 -jit -var stim_period 2000.00 -var T 298. -downsampleRate 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.05 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p21_serca1p00_T298_freq0p5_dc.pickle &

#HIP (SERCA=50%)
python daisychain.py -dt 0.1 -jit -var stim_period 2000.00 -var T 298. -downsampleRate 10. -iters 30 -T 10000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.05 -var V_max_Jpump 0.00351 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p21_serca0p50_freq0p5_dc.pickle &

#2.0 Hz
#Amylin
#SERCA=100%
python daisychain.py -dt 0.1 -jit -var stim_period 500.00 -var T 298. -downsampleRate 10. -iters 30 -T 10000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.05 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p21_serca1p00_T298_freq2p0_dc.pickle &

#HIP (SERCA=50%)
python daisychain.py -dt 0.1 -jit -var stim_period 500.00 -var T 298. -downsampleRate 10. -iters 30 -T 10000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.05 -var V_max_Jpump 0.00351 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p21_serca0p50_T298_freq2p0_dc.pickle &

#0.5 Hz
#WT
#SERCA=100%
#python daisychain.py -dt 0.1 -jit -var stim_period 2000.00 -var T 298. -downsampleRate 10. -iters 30 -T 10000 -odeName shannon_2004_mouse.ode -name mouse_leak1p00_nka1p00_serca1p00_T298_freq0p5_dc.pickle &

#2.0 Hz
#WT
#SERCA=100%
python daisychain.py -dt 0.1 -jit -var stim_period 500.00 -var T 298. -downsampleRate 10. -iters 30 -T 10000 -odeName shannon_2004_mouse.ode -name mouse_leak1p00_nka1p00_serca1p00_T298_freq2p0_dc.pickle &

