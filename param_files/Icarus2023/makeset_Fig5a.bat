krn="P300"
echo "$krn""_lon0">listnames.txt
echo "$krn""_lon20">>listnames.txt
echo "$krn""_lon40">>listnames.txt
echo "$krn""_lon60">>listnames.txt
echo "$krn""_lon80">>listnames.txt
echo "$krn""_lon90">>listnames.txt
xargs mkdir <listnames.txt

for folder in ./P*/
do
	cd $folder
	cp -r ../src/* ./
	rm ./param_Fig5b.in
	mv ./param_Fig5a.in ./param.in
	tmpvar=$(cut -d "n" -f2 <<< "$folder")
	tmpvar=$(head -c 2 <<< "$tmpvar")
	echo "folder $folder, cap_lon $tmpvar"
	sed -i "s/cap_lon     =0./cap_lon     =$tmpvar./g" ./param.in
	cd ../
done

for folder in ./P*/
do
	tmpvar=${folder%/*}
	tmpvar=${tmpvar##*/}
    mkdir "$tmpvar"_fb
    cp -r ./$tmpvar/* ./"$tmpvar"_fb/
    sed -i 's/fosslit     =.false./fosslit     =.true./g' ./"$tmpvar"_fb/param.in
done
rm listnames.txt
