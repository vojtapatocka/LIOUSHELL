proc=$(grep loadstate ./param.in | awk '{print $2}')
echo "Postprocessing script, you may provide id (default 11), compdir (default pure), np (to avoid the polar plot)"

if [ -f "./figs/balance_pr0_pure.png" ]
then
    proc="=1"
fi
if [ -z "$1" ]
then
    id=10
    tides=$(grep tides ./param.in | awk '{print $2}')
    if [ "$tides" == "=.true." ]
    then
        id=11
    fi
else
    id=$1
fi
if [ -z "$2" ]
then
    compdir='pure'
else
    sid=$(echo "$2" | sed 's=../==g')
    compdir=$sid
fi
echo "Running gnuplot scripts with id=$id and compdir=$compdir"

mkdir figs

if [ -f "./urad_ini.tail" ]
then
    Erot0=$(tail -n 1 "./urad_ini.tail" | awk '{print $4}')
    Etid0=$(tail -n 1 "./urad_ini.tail" | awk '{print $8}')
    Egrav0=$(tail -n 1 "./urad_ini.tail" | awk '{print $5}')
    Eel0=$(tail -n 1 "./urad_ini.tail" | awk '{print $6}')
else
    Erot0=0; Egrav0=0; Eel0=0; Etid0=0;
fi

cd gnuplot
gnuplot<profiles.g
if [ "$proc" == "=1" ]
then
    echo "postprocessing: TPW (pr=1)"
    if [ "$id" == "p" ]
    then
        echo "Running only polar.py, compdir=$2"
        python polar.py "compdir='$2'"
    else
        echo "full TPW postprocessing: polar.py and gnuplot scripts (spada, solution, terms)"
gnuplot <<- EOF
    pr="1"
    Erot0="${Erot0}"
    Etid0="${Etid0}"
    Egrav0="${Egrav0}"
    Eel0="${Eel0}"
    id="${id}"
    compdir="${compdir}"
    l 'spada.g'
    l 'solution.g'
    l 'terms.g'
EOF
        if [ "$3" != "np" ] ; then python polar.py "compdir='$compdir'"; fi
    fi
else
    echo "postprocessing: RELAXATION (pr=0)"
gnuplot <<- EOF
    pr="0"
    Erot0="${Erot0}"
    Etid0="${Etid0}"
    Egrav0="${Egrav0}"
    Eel0="${Eel0}"
    compdir="${compdir}"
    l 'solution.g'
    l 'terms.g'
EOF
fi

cd ../
