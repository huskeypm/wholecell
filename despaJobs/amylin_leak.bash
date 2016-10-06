#Try to replicate Despa findings.
#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -dSr 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p00_freq1p0_dc.pickle &

#Increase NKA
#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -dSr 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.0 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p20_freq1p0_dc.pickle &

#Amylin
#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -downsampleRate 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.05 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p21_hip1p00_freq1p0_dc.pickle &

#0.5 Hz
#python daisychain.py -dt 0.1 -jit -var stim_period 2000.00 -downsampleRate 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.05 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p21_hip0p50_freq0p50_dc.pickle &

#2.0 Hz
#python daisychain.py -dt 0.1 -jit -var stim_period 500.00 -downsampleRate 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.05 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p21_hip0p50_freq2p00_dc.pickle &

#1.0 Hz, change NKA function
#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -dSr 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.1 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p22_freq1p0_dc.pickle &

#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -dSr 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.125 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p225_freq1p0_dc.pickle &

#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -dSr 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.15 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p23_freq1p0_dc.pickle &

#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -dSr 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.25 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p25_freq1p0_dc.pickle &

#WT 5min
#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -downsampleRate 10. -iters 15 -T 20000 -odeName shannon_2004_mouse.ode -name mouse_leak1p000_nka1p00_hip1p00_freq1p0_dc.pickle &

#HIP
#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -downsampleRate 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.05 -var V_max_Jpump 0.00351 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p21_hip0p50_freq1p0_dc.pickle &

#0.5 Hz
#python daisychain.py -dt 0.1 -jit -var stim_period 2000.00 -downsampleRate 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.05 -var V_max_Jpump 0.00351 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p21_hip0p50_freq0p50_dc.pickle &

#2.0 Hz
#python daisychain.py -dt 0.1 -jit -var stim_period 500.00 -downsampleRate 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -var I_NaK_max 6.05 -var V_max_Jpump 0.00351 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p21_hip0p50_freq2p00_dc.pickle &

#WT 
#0.5HZ
#python daisychain.py -dt 0.1 -jit -var stim_period 2000.00 -downsampleRate 10. -iters 15 -T 20000 -odeName shannon_2004_mouse.ode -name mouse_leak1p000_nka1p00_hip1p00_freq0p50_dc.pickle &

#2.0 Hz
#python daisychain.py -dt 0.1 -jit -var stim_period 500.00 -downsampleRate 10. -iters 15 -T 20000 -odeName shannon_2004_mouse.ode -name mouse_leak1p000_nka1p00_hip1p00_freq2p00_dc.pickle &

##TEST
#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 --downsampleRate 10. -iters 6 -T 2000 -odeName shannon_2004_mouse.ode -name temp10_dc.pickle &

#NKA cases
#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -downsampleRate 10. -iters 15 -T 20000 -var I_NaK_max 6.05 -odeName shannon_2004_mouse.ode -name mouse_leak1p00_nka1p21_freq1p0_dc.pickle &

#Leak cases
#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -downsampleRate 10. -iters 15 -T 20000 -var G_CaBk 0.0010743075 -var G_NaBk 0.001905225 -odeName shannon_2004_mouse.ode -name mouse_leak1p425_nka1p00_hip1p00_freq1p0_dc.pickle &

#Test cases

#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -dSr 10. -iters 2 -T 2000 -var I_NaK_max 6.05 -odeName shannon_2004_mouse.ode -name test_daisychain.pickle 

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.001319325 -T 30000 -odeName shannon_2004_mouse.ode -name test_runShannon.pickle &

#LCC
python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -downsampleRate 10. -iters 15 -T 20000 -var PCa 0.00111375 -odeName shannon_2004_mouse.ode -name mouse_leak1p00x_lcc1p25_freq1p0_dc.pickle &
