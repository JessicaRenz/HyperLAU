g++ HyperLAU.cpp -o HyperLAU -larmadillo

set -- 

#run examples from the data folder
#First toy example without bootstrap, random seed 1, model F
./HyperLAU data/first_toyexample.txt 3 first_toyexample_mF 0 1 -1 1.001 > first_toyexample_mF.tmp &
set -- "$@" $!

#First toy example without bootstrap, random seed 1, model 1
./HyperLAU data/first_toyexample.txt 3 first_toyexample_m1 0 1 1 1.001 > first_toyexample_m1.tmp &
set -- "$@" $!

#First toy example without bootstrap, random seed 1, model 2
./HyperLAU data/first_toyexample.txt 3 first_toyexample_m2 0 1 2 1.001 > first_toyexample_m2.tmp &
set -- "$@" $!

#Second toy example with inserted uncertainty markers in different features, all with model F
./HyperLAU data/second_toyexample/full.txt 6 full_mF 0 1 -1 1.001 > full_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature1_40.txt 6 feature1_40_mF 0 1 -1 1.001 > feature1_40_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature2_40.txt 6 feature2_40_mF 0 1 -1 1.001 > feature2_40_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature3_40.txt 6 feature3_40_mF 0 1 -1 1.001 > feature3_40_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature4_40.txt 6 feature4_40_mF 0 1 -1 1.001 > feature4_40_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature5_40.txt 6 feature5_40_mF 0 1 -1 1.001 > feature5_40_mF.tmp &
set -- "$@" $!
./HyperLAU data/second_toyexample/feature6_40.txt 6 feature6_40_mF 0 1 -1 1.001 > feature6_40_mF.tmp &
set -- "$@" $!

#the tuberculosis datasets, both with model F
./HyperLAU data/tb/tb_data_10.txt 10 tb_data_10_mF 0 1 -1 1.001 > tb_data_10_mF.tmp &
set -- "$@" $!
./HyperLAU data/tb/tb_data_10_qm50.txt 10 tb_data_10_qm50_mF 0 1 -1 1.001 > tb_data_10_qm50_mF.tmp &
set -- "$@" $!

echo PIDS: $@
for job in $@
do
    wait $job
    echo $job finished
done
