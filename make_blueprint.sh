for i in $1
do

for j in $number_of_individuals_for_each_location
do

echo "#input setting"
echo popid: ${i}_allsites_vartest
echo nseq: $((j*2))
echo L: 4000000
echo whether_folded: true
echo SFS: `cat ${i}_sfs`
echo dimFactor: 5000 "# dimension factor for numerical optimization"
echo pct_training: 0.67 "# percentage of sites for training"
echo nrand: $((j/2-1/2))        $((j-1))      $((j*3/2-3/2))      $((j*4/2-4/2)) "# number of random break points for each try (separated by white space)"
echo project_dir: ${i}_allsites_vartest "# project directory"
echo stairway_plot_dir: stairway_plot_es "# directory to the stairway plot files"
echo ninput: 200 "# number of input files to be created for each estimation"
echo "#output setting"
echo mu: 1e-7 "# assumed mutation rate per site per generation"
echo year_per_generation: 5 "# assumed generation time (in years)"
echo "#plot setting"
echo plot_title: ${i}_allsites_vartest "# title of the plot"
echo xrange: 0,0 "# Time (1k year) range; format: xmin,xmax; "0,0" for default"
echo yrange: 0,0 "# Ne (1k individual) range; format: xmin,xmax; "0,0" for default"
echo xspacing: 2 "# X axis spacing"
echo yspacing: 2 "# Y axis spacing"
echo fontsize: 14 "# Font size"


done

done
