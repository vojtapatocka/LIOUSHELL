ln="lon4"
echo "N030_$ln">listnames.txt
echo "N100_$ln">>listnames.txt
echo "N200_$ln">>listnames.txt
echo "N200_lon0">>listnames.txt
echo "N300_$ln">>listnames.txt
xargs mkdir <listnames.txt

for folder in ./N*/
do
	cd $folder
	cp -r ../src/* ./
	rm ./param_Fig5a.in
	mv ./param_Fig5b.in ./param.in
	tmpvar=$(cut -d "N" -f2 <<< "$folder")
	tmpvar=$(head -c 3 <<< "$tmpvar")
	echo "folder $folder, hice $tmpvar"
	sed -i "s/hice        =200/hice        =$tmpvar/g" ./param.in
	cd ../
done

for folder in ./N*/
do
	tmpvar=${folder%/*}
	tmpvar=${tmpvar##*/}
    mkdir "$tmpvar"_fb
    cp -r ./$tmpvar/* ./"$tmpvar"_fb/
    sed -i 's/fosslit     =.false./fosslit     =.true./g' ./"$tmpvar"_fb/param.in
done
sed -i 's/cap_lon     =0./cap_lon     =4./g' ./*lon4*/param.in
