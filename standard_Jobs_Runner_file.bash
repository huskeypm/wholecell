### control, Amy, and HIP over freq 1.0 Hz

#### making control data ####

python daisychain.py -dt 0.10 -jit -stim_period 1000 -T 10000 -iters 30 -fileOutputDirectory /home/AD/bdst227/ipython/ipython-notebooks/Despa/wholecell/despaJobs/ranJobs/ -finalOutputDirectory /net/share/bdst227/Despa/Despa_Simulations_Data/ -downsampleRate 10 -state Na_SL 12.00 -state Na_jct1 12.00 -state Nai 12.00 -var G_CaBk 0.0007539 -var G_NaBk 0.001337 -var I_NaK_max 3.85 -var T 310.00 -var V_max_Jpump 0.009977826 -odeName shannon_2004_rat.ode -name /home/AD/bdst227/ipython/ipython-notebooks/Despa/wholecell/despaJobs/ranJobs/rat_BASELINE_freq1p0Hz &

#### making Amy data ####

python daisychain.py -dt 0.10 -jit -stim_period 1000 -T 10000 -iters 30 -fileOutputDirectory /home/AD/bdst227/ipython/ipython-notebooks/Despa/wholecell/despaJobs/ranJobs/ -finalOutputDirectory /net/share/bdst227/Despa/Despa_Simulations_Data/ -downsampleRate 10 -state Na_SL 12.00 -state Na_jct1 12.00 -state Nai 12.00 -var G_CaBk 0.001267 -var G_NaBk 0.001337 -var I_NaK_max 4.4 -var T 310.00 -var V_max_Jpump 0.009977826 -odeName shannon_2004_rat.ode -name /home/AD/bdst227/ipython/ipython-notebooks/Despa/wholecell/despaJobs/ranJobs/rat_AMY_freq1p0Hz &

#### making HIP data ####

python daisychain.py -dt 0.10 -jit -stim_period 1000 -T 10000 -iters 30 -fileOutputDirectory /home/AD/bdst227/ipython/ipython-notebooks/Despa/wholecell/despaJobs/ranJobs/ -finalOutputDirectory /net/share/bdst227/Despa/Despa_Simulations_Data/ -downsampleRate 10 -state Na_SL 12.00 -state Na_jct1 12.00 -state Nai 12.00 -var G_CaBk 0.0009424 -var G_NaBk 0.001337 -var I_NaK_max 4.02 -var T 310.00 -var V_max_Jpump 0.006156 -odeName shannon_2004_rat.ode -name /home/AD/bdst227/ipython/ipython-notebooks/Despa/wholecell/despaJobs/ranJobs/rat_HIP_freq1p0Hz &

