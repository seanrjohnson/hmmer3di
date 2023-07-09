NEWDIR="hmmer3Di"
prefix="3Di_"

mkdir -p $NEWDIR
for f in hmmalign hmmbuild hmmpress hmmsearch hmmscan phmmer
do
    cp src/${f} ${NEWDIR}/${prefix}${f}
done

