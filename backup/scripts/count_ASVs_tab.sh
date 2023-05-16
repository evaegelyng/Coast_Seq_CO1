
gwf status tab* | grep "shouldrun" > results/tab_shouldrun.txt
cut -d " " -f1 results/tab_shouldrun.txt | sed "s/tab_//g" > results/motus_shouldrun.txt
file="results/motus_shouldrun.txt"
rm_force results/asv_count.txt
for name in $(cat $file)
do
    count=$(wc -l tmp/motus/$name | cut -d " " -f1)
    echo $name $count >> results/asv_count.txt
done