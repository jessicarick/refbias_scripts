#!/bin/sh

phy=$1
base=`echo $phy | sed 's/\.phy//g'`

nind=$(head -n 1 $phy | cut -f 1 -d' ')
nsites=$(head -n 1 $phy | cut -f 2 -d' ')

echo "working with phylip with $nind individuals and $nsites sites"

rm -f tmp.rmv

# count number of N per line
i=2
n=0
tail -n +2 $phy | cut -f 2 | awk -F'N' '{print NF-1}' | while read nN; do
	if [[ nN -ge nsites ]]; then
		if [[ n -eq 0 ]]; then
			echo "${i}d" > tmp.rmv
		else
			echo ";${i}d" >> tmp.rmv
		fi

		i=$((i+1))
		n=$((n+1))
	else
		i=$((i+1))
	fi
done 

ndelete=`cat tmp.rmv | wc -l`

if [[ ndelete -ge 1 ]]; then
	delete=`cat tmp.rmv | tr -d '\n'`
        ndelete=`cat tmp.rmv | wc -l`
        new_nind=$((nind - ndelete))

        echo "removing $ndelete individuals that are completely missing"
        echo "$new_nind $nsites" > ${base}.reduced.phy
        sed -e $delete $phy | tail -n +2 >> ${base}.reduced.phy
else
	echo "no individuals to remove"
	cp $phy ${base}.reduced.phy
fi
