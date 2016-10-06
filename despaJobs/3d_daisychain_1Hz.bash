#Amylin
python daisychain.py -dt 0.1 -jit -var G_CaBk 0.000942375 -var G_NaBk 0.00167125 -var stim_period 1000.00 -dSr 10. -iters 15 -T 20000 -odeName shannon_2004_mouse.ode -name mouse_leak1p25x_SERCAuptake1p00x_freq1p0_dc.pickle &

#HIP
python daisychain.py -dt 0.1 -jit -var G_CaBk 0.000942375 -var G_NaBk 0.00167125 -var V_max_Jpump 0.003510 -var stim_period 1000.00 -dSr 10. -iters 15 -T 20000 -odeName shannon_2004_mouse.ode -name mouse_leak1p25x_SERCAuptake0p50x_freq1p0_dc.pickle &

#UCD
#python daisychain.py -dt 0.1 -jit -var V_max_Jpump 0.003510 -var stim_period 1000.00 -dSr 10. -iters 5 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p00x_SERCAuptake0p50x_freq1p0_dc.pickle &

#Increase leak

#python daisychain.py -dt 0.1 -jit -var stim_period 1000.00 -dSr 10. -iters 5 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p00x_freq1p0_dc.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -T 60000 -odeName shannon_2004_mouse.ode -name tmp.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00098007 -var G_NaBk 0.0017381 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00105546 -var G_NaBk 0.0018718 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00067851 -var G_NaBk 0.0012033 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00060312 -var G_NaBk 0.0010696 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00052773 -var G_NaBk 0.0009359 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00045234 -var G_NaBk 0.0008022 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00037695 -var G_NaBk 0.0006685 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p50x_freq1p0_60.pickle &

#Decrease SERCA uptake

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var V_max_Jpump 0.003510 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_SERCAuptake0p50x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var V_max_Jpump 0.004212 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_SERCAuptake0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var V_max_Jpump 0.004914 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_SERCAuptake0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var V_max_Jpump 0.005616 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_SERCAuptake0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var V_max_Jpump 0.006318 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_SERCAuptake0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var V_max_Jpump 0.007722 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_SERCAuptake1p10x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var V_max_Jpump 0.008424 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_SERCAuptake1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var V_max_Jpump 0.009126 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_SERCAuptake1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var V_max_Jpump 0.009828 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_SERCAuptake1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var V_max_Jpump 0.010530 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_SERCAuptake1p50x_freq1p0_60.pickle &

# Leak=50%, Change SERCA
#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00037695 -var G_NaBk 0.0006685 -var V_max_Jpump 0.003510 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p50x_SERCAuptake0p50x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00037695 -var G_NaBk 0.0006685 -var V_max_Jpump 0.004212 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p50x_SERCAuptake0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00037695 -var G_NaBk 0.0006685 -var V_max_Jpump 0.004914 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p50x_SERCAuptake0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00037695 -var G_NaBk 0.0006685 -var V_max_Jpump 0.005616 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p50x_SERCAuptake0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00037695 -var G_NaBk 0.0006685 -var V_max_Jpump 0.006318 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p50x_SERCAuptake0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00037695 -var G_NaBk 0.0006685 -var V_max_Jpump 0.007722 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p50x_SERCAuptake1p10x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00037695 -var G_NaBk 0.0006685 -var V_max_Jpump 0.008424 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p50x_SERCAuptake1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00037695 -var G_NaBk 0.0006685 -var V_max_Jpump 0.009126 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p50x_SERCAuptake1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00037695 -var G_NaBk 0.0006685 -var V_max_Jpump 0.009828 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p50x_SERCAuptake1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00037695 -var G_NaBk 0.0006685 -var V_max_Jpump 0.010530 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p50x_SERCAuptake1p50x_freq1p0_60.pickle &

# Leak=60%, Change SERCA
#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00045234 -var G_NaBk 0.0008022 -var V_max_Jpump 0.003510 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p60x_SERCAuptake0p50x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00045234 -var G_NaBk 0.0008022 -var V_max_Jpump 0.004212 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p60x_SERCAuptake0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00045234 -var G_NaBk 0.0008022 -var V_max_Jpump 0.004914 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p60x_SERCAuptake0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00045234 -var G_NaBk 0.0008022 -var V_max_Jpump 0.005616 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p60x_SERCAuptake0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00045234 -var G_NaBk 0.0008022 -var V_max_Jpump 0.006318 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p60x_SERCAuptake0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00045234 -var G_NaBk 0.0008022 -var V_max_Jpump 0.007722 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p60x_SERCAuptake1p10x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00045234 -var G_NaBk 0.0008022 -var V_max_Jpump 0.008424 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p60x_SERCAuptake1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00045234 -var G_NaBk 0.0008022 -var V_max_Jpump 0.009126 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p60x_SERCAuptake1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00045234 -var G_NaBk 0.0008022 -var V_max_Jpump 0.009828 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p60x_SERCAuptake1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00045234 -var G_NaBk 0.0008022 -var V_max_Jpump 0.010530 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p60x_SERCAuptake1p50x_freq1p0_60.pickle &

# Leak=70%, Change SERCA
#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00052773 -var G_NaBk 0.0009359 -var V_max_Jpump 0.003510 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p70x_SERCAuptake0p50x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00052773 -var G_NaBk 0.0009359 -var V_max_Jpump 0.004212 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p70x_SERCAuptake0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00052773 -var G_NaBk 0.0009359 -var V_max_Jpump 0.004914 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p70x_SERCAuptake0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00052773 -var G_NaBk 0.0009359 -var V_max_Jpump 0.005616 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p70x_SERCAuptake0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00052773 -var G_NaBk 0.0009359 -var V_max_Jpump 0.006318 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p70x_SERCAuptake0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00052773 -var G_NaBk 0.0009359 -var V_max_Jpump 0.007722 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p70x_SERCAuptake1p10x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00052773 -var G_NaBk 0.0009359 -var V_max_Jpump 0.008424 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p70x_SERCAuptake1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00052773 -var G_NaBk 0.0009359 -var V_max_Jpump 0.009126 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p70x_SERCAuptake1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00052773 -var G_NaBk 0.0009359 -var V_max_Jpump 0.009828 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p70x_SERCAuptake1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00052773 -var G_NaBk 0.0009359 -var V_max_Jpump 0.010530 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p70x_SERCAuptake1p50x_freq1p0_60.pickle &

# Leak=80%, Change SERCA
#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00060312 -var G_NaBk 0.0010696 -var V_max_Jpump 0.003510 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p80x_SERCAuptake0p50x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00060312 -var G_NaBk 0.0010696 -var V_max_Jpump 0.004212 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p80x_SERCAuptake0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00060312 -var G_NaBk 0.0010696 -var V_max_Jpump 0.004914 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p80x_SERCAuptake0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00060312 -var G_NaBk 0.0010696 -var V_max_Jpump 0.005616 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p80x_SERCAuptake0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00060312 -var G_NaBk 0.0010696 -var V_max_Jpump 0.006318 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p80x_SERCAuptake0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00060312 -var G_NaBk 0.0010696 -var V_max_Jpump 0.007722 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p80x_SERCAuptake1p10x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00060312 -var G_NaBk 0.0010696 -var V_max_Jpump 0.008424 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p80x_SERCAuptake1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00060312 -var G_NaBk 0.0010696 -var V_max_Jpump 0.009126 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p80x_SERCAuptake1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00060312 -var G_NaBk 0.0010696 -var V_max_Jpump 0.009828 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p80x_SERCAuptake1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00060312 -var G_NaBk 0.0010696 -var V_max_Jpump 0.010530 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p80x_SERCAuptake1p50x_freq1p0_60.pickle &

# Leak=90%, Change SERCA
#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00067851 -var G_NaBk 0.0012033 -var V_max_Jpump 0.003510 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p90x_SERCAuptake0p50x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00067851 -var G_NaBk 0.0012033 -var V_max_Jpump 0.004212 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p90x_SERCAuptake0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00067851 -var G_NaBk 0.0012033 -var V_max_Jpump 0.004914 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p90x_SERCAuptake0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00067851 -var G_NaBk 0.0012033 -var V_max_Jpump 0.005616 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p90x_SERCAuptake0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00067851 -var G_NaBk 0.0012033 -var V_max_Jpump 0.006318 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p90x_SERCAuptake0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00067851 -var G_NaBk 0.0012033 -var V_max_Jpump 0.007722 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p90x_SERCAuptake1p10x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00067851 -var G_NaBk 0.0012033 -var V_max_Jpump 0.008424 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p90x_SERCAuptake1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00067851 -var G_NaBk 0.0012033 -var V_max_Jpump 0.009126 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p90x_SERCAuptake1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00067851 -var G_NaBk 0.0012033 -var V_max_Jpump 0.009828 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p90x_SERCAuptake1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00067851 -var G_NaBk 0.0012033 -var V_max_Jpump 0.010530 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak0p90x_SERCAuptake1p50x_freq1p0_60.pickle &

# Leak=110%, Change SERCA
#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00082929 -var G_NaBk 0.0014707 -var V_max_Jpump 0.003510 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p10x_SERCAuptake0p50x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00082929 -var G_NaBk 0.0014707 -var V_max_Jpump 0.004212 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p10x_SERCAuptake0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00082929 -var G_NaBk 0.0014707 -var V_max_Jpump 0.004914 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p10x_SERCAuptake0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00082929 -var G_NaBk 0.0014707 -var V_max_Jpump 0.005616 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p10x_SERCAuptake0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00082929 -var G_NaBk 0.0014707 -var V_max_Jpump 0.006318 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p10x_SERCAuptake0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00082929 -var G_NaBk 0.0014707 -var V_max_Jpump 0.007722 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p10x_SERCAuptake1p10x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00082929 -var G_NaBk 0.0014707 -var V_max_Jpump 0.008424 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p10x_SERCAuptake1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00082929 -var G_NaBk 0.0014707 -var V_max_Jpump 0.009126 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p10x_SERCAuptake1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00082929 -var G_NaBk 0.0014707 -var V_max_Jpump 0.009828 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p10x_SERCAuptake1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00082929 -var G_NaBk 0.0014707 -var V_max_Jpump 0.010530 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p10x_SERCAuptake1p50x_freq1p0_60.pickle &

# Leak=120%, Change SERCA
#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -var V_max_Jpump 0.003510 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p20x_SERCAuptake0p50x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -var V_max_Jpump 0.004212 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p20x_SERCAuptake0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -var V_max_Jpump 0.004914 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p20x_SERCAuptake0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -var V_max_Jpump 0.005616 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p20x_SERCAuptake0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -var V_max_Jpump 0.006318 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p20x_SERCAuptake0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -var V_max_Jpump 0.007722 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p20x_SERCAuptake1p10x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -var V_max_Jpump 0.008424 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p20x_SERCAuptake1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -var V_max_Jpump 0.009126 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p20x_SERCAuptake1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -var V_max_Jpump 0.009828 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p20x_SERCAuptake1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00090468 -var G_NaBk 0.0016044 -var V_max_Jpump 0.010530 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p20x_SERCAuptake1p50x_freq1p0_60.pickle &

# Leak=130%, Change SERCA
#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00098007 -var G_NaBk 0.0017381 -var V_max_Jpump 0.003510 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p30x_SERCAuptake0p50x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00098007 -var G_NaBk 0.0017381 -var V_max_Jpump 0.004212 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p30x_SERCAuptake0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00098007 -var G_NaBk 0.0017381 -var V_max_Jpump 0.004914 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p30x_SERCAuptake0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00098007 -var G_NaBk 0.0017381 -var V_max_Jpump 0.005616 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p30x_SERCAuptake0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00098007 -var G_NaBk 0.0017381 -var V_max_Jpump 0.006318 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p30x_SERCAuptake0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00098007 -var G_NaBk 0.0017381 -var V_max_Jpump 0.007722 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p30x_SERCAuptake1p10x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00098007 -var G_NaBk 0.0017381 -var V_max_Jpump 0.008424 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p30x_SERCAuptake1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00098007 -var G_NaBk 0.0017381 -var V_max_Jpump 0.009126 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p30x_SERCAuptake1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00098007 -var G_NaBk 0.0017381 -var V_max_Jpump 0.009828 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p30x_SERCAuptake1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00098007 -var G_NaBk 0.0017381 -var V_max_Jpump 0.010530 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p30x_SERCAuptake1p50x_freq1p0_60.pickle &

# Leak=140%, Change SERCA
#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00105546 -var G_NaBk 0.0018718 -var V_max_Jpump 0.003510 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p40x_SERCAuptake0p50x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00105546 -var G_NaBk 0.0018718 -var V_max_Jpump 0.004212 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p40x_SERCAuptake0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00105546 -var G_NaBk 0.0018718 -var V_max_Jpump 0.004914 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p40x_SERCAuptake0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00105546 -var G_NaBk 0.0018718 -var V_max_Jpump 0.005616 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p40x_SERCAuptake0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00105546 -var G_NaBk 0.0018718 -var V_max_Jpump 0.006318 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p40x_SERCAuptake0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00105546 -var G_NaBk 0.0018718 -var V_max_Jpump 0.007722 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p40x_SERCAuptake1p10x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00105546 -var G_NaBk 0.0018718 -var V_max_Jpump 0.008424 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p40x_SERCAuptake1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00105546 -var G_NaBk 0.0018718 -var V_max_Jpump 0.009126 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p40x_SERCAuptake1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00105546 -var G_NaBk 0.0018718 -var V_max_Jpump 0.009828 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p40x_SERCAuptake1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00105546 -var G_NaBk 0.0018718 -var V_max_Jpump 0.010530 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p40x_SERCAuptake1p50x_freq1p0_60.pickle &

# Leak=150%, Change SERCA
#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00113085 -var G_NaBk 0.0020055 -var V_max_Jpump 0.003510 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p50x_SERCAuptake0p50x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00113085 -var G_NaBk 0.0020055 -var V_max_Jpump 0.004212 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p50x_SERCAuptake0p60x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00113085 -var G_NaBk 0.0020055 -var V_max_Jpump 0.004914 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p50x_SERCAuptake0p70x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00113085 -var G_NaBk 0.0020055 -var V_max_Jpump 0.005616 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p50x_SERCAuptake0p80x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00113085 -var G_NaBk 0.0020055 -var V_max_Jpump 0.006318 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p50x_SERCAuptake0p90x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00113085 -var G_NaBk 0.0020055 -var V_max_Jpump 0.007722 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p50x_SERCAuptake1p10x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00113085 -var G_NaBk 0.0020055 -var V_max_Jpump 0.008424 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p50x_SERCAuptake1p20x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00113085 -var G_NaBk 0.0020055 -var V_max_Jpump 0.009126 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p50x_SERCAuptake1p30x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00113085 -var G_NaBk 0.0020055 -var V_max_Jpump 0.009828 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p50x_SERCAuptake1p40x_freq1p0_60.pickle &

#python runShannonTest.py -dt 0.1 -stim 1000.00 -jit -var G_CaBk 0.00113085 -var G_NaBk 0.0020055 -var V_max_Jpump 0.010530 -T 60000 -odeName shannon_2004_mouse.ode -name mouse_leak1p50x_SERCAuptake1p50x_freq1p0_60.pickle &
