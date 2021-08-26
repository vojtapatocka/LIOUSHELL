gnuplot<solution.g
gnuplot<terms.g
gnuplot<profiles.g
gnuplot<love_number.g

sed -i 's/pr=0/pr=1/g' ./solution.g
sed -i 's/pr=0/pr=1/g' ./terms.g

gnuplot<solution.g
gnuplot<terms.g

sed -i 's/pr=1/pr=0/g' ./solution.g
sed -i 's/pr=1/pr=0/g' ./terms.g
