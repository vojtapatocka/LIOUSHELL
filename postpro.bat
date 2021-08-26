proc=$(grep loadstate ./param.in | awk '{print $2}')

if [ -f "./figs/balance_pr0_pure.png" ]
then
    proc="=1"
fi
if [ -z "$1" ]
then
    id=10
else
    id=$1
fi
if [ -z "$2" ]
then
    compdir='pure'
    echo "Not running a comparison, disregard warnings for missing ../../pure files"
else
    sid=$(echo "$2" | sed 's=../==g')
    compdir=$sid
fi
echo "Running gnuplot scripts with id=$id and compdir=$compdir"

mkdir figs

if [ -f "../urad_ini.tail" ]
then
    Erot0=$(tail -n 1 "../urad_ini.tail" | awk '{print $4}')
    Egrav0=$(tail -n 1 "../urad_ini.tail" | awk '{print $5}')
    Eel0=$(tail -n 1 "../urad_ini.tail" | awk '{print $6}')
else
    Erot0=0
    Egrav0=0
    Eel0=0
fi

cd gnuplot
gnuplot<profiles.g
if [ "$proc" == "=1" ]
then
gnuplot <<- EOF
    pr="1"
    Erot0="${Erot0}"
    Egrav0="${Egrav0}"
    Eel0="${Eel0}"
    id="${id}"
    compdir="${compdir}"
    l 'spada.g'
    l 'solution.g'
EOF
else
gnuplot <<- EOF
    pr="0"
    Erot0="${Erot0}"
    Egrav0="${Egrav0}"
    Eel0="${Eel0}"
    compdir="${compdir}"
    l 'solution.g'
EOF
fi

cd ../

