echo "N100_hfP">listnames.txt
echo "N200_hfP">>listnames.txt
echo "N300_hfP">>listnames.txt
echo "N300_hf001">>listnames.txt
echo "N300_hf1">>listnames.txt
echo "N300_hfinf">>listnames.txt
xargs mkdir <listnames.txt

for folder in ./N*/
do
	cd $folder
	cp -r ../src/* ./
	tmpvar=$(cut -d "N" -f2 <<< "$folder")
	tmpvar=$(head -c 3 <<< "$tmpvar")
	echo "folder $folder, hice $tmpvar"
	sed -i "s/hice        =300/hice        =$tmpvar/g" ./param.in
	cd ../
done
sed -i 's/0.121161192/0.01/g' ./*hf001/param.in
sed -i 's/0.121161192/1.0/g' ./*hf1/param.in
sed -i 's/0.121161192/1.0e10/g' ./*hfinf/param.in

for folder in ./N*/
do
	tmpvar=${folder%/*}
	tmpvar=${tmpvar##*/}
	pvar=$(sed 's/N/P/g' <<< $tmpvar)
	echo "making $pvar directory"
    mkdir "$pvar"
    cp -r ./$tmpvar/* ./$pvar/
    sed -i 's/cap_col     =89/cap_col     =1/g' ./$pvar/param.in    
    sed -i 's/rhoice      =-1000/rhoice      =1000/g' ./$pvar/param.in
done

for folder in ./*hf*/
do
	tmpvar=${folder%/*}
	tmpvar=${tmpvar##*/}
    mkdir "$tmpvar"_fb
    cp -r ./$tmpvar/* ./"$tmpvar"_fb/
    sed -i 's/fosslit     =.false./fosslit     =.true./g' ./"$tmpvar"_fb/param.in
done
rm listnames.txt
