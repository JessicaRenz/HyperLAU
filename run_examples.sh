g++ HyperLAU.cpp -o HyperLAU -larmadillo

set -- 

#run examples from the data folder
#First toy example without bootstrap, random seed 1, model F
./HyperLAU data/first_toyexample.txt first_toyexample_mF > first_toyexample_mF.tmp &
set -- "$@" $!

#First toy example without bootstrap, random seed 1, model 1
./HyperLAU data/first_toyexample.txt first_toyexample_m1 --model 1 > first_toyexample_m1.tmp &
set -- "$@" $!

#First toy example without bootstrap, random seed 1, model 2
./HyperLAU data/first_toyexample.txt first_toyexample_m2 --model 2 > first_toyexample_m2.tmp &
set -- "$@" $!

#Second toy example with inserted uncertainty markers in different features, all with model F
./HyperLAU data/second_toyexample/full.txt full_mF > full_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature1_40.txt feature1_40_mF > feature1_40_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature2_40.txt feature2_40_mF > feature2_40_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature3_40.txt feature3_40_mF > feature3_40_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature4_40.txt feature4_40_mF > feature4_40_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature5_40.txt feature5_40_mF > feature5_40_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature6_40.txt feature6_40_mF > feature6_40_mF.tmp &
set -- "$@" $!

#the tuberculosis datasets, both with model F
./HyperLAU data/tb/tb_data_10.txt tb_data_10_mF > tb_data_10_mF.tmp &
set -- "$@" $!
./HyperLAU data/tb/tb_data_10_qm50.txt tb_data_10_qm50_mF > tb_data_10_qm50_mF.tmp &
set -- "$@" $!

echo PIDS: $@
for job in $@
do
    wait $job
    echo $job finished
done
