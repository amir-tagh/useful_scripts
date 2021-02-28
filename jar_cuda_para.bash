#!/bin/bash -l

##load the modules command
. /etc/profile.d/resif.sh
. /etc/profile.d/nvidia.sh


#module load OpenMPI/1.7.3-gcccuda-2.6.10 
module load numlib/imkl
#module load CUDA 

module load toolchain/iimpi
module load toolchain/ictce
module load mpi/impi

stem=$1
card=$2


##build restraint and infiles
../make_rst.bash $stem.pdb 

##meant to be run whole-node


get_free_card () {
    cards_present=$(nvidia-smi -L | awk '{print $2}' | sed 's/://')
    cards_used=$(nvidia-smi  | grep pmemd | awk '{print $2}')
    node=$(cat $OAR_NODEFILE | head -n 1)
   use_c=-1
  for c in $cards_present
  do   
     ok=1
     for cc in $cards_used
     do    
           if [ "$cc" -eq "$c" ]
           then
                   ok=0
           fi
     done
     if [ "$ok" -eq "1" ]
     then
           use_c=$c
           break
     else
           echo "card $c is busy"
     fi
  done
  if [ "$use_c" -eq "-1" ]
  then
    echo "did not find non-busy GPU on node $node"
    echo "consider a resub asking for access to an entire node instead of 1 core."
    nvidia-smi
    exit 1
  else
    echo "running on non-busy GPU, id: $use_c"
    export CUDA_VISIBLE_DEVICES=$use_c
  fi

  ##return total number available... just set a variale, its global.
  count_cards_present=$(echo $cards_present | wc | awk '{print $2}')
}

##call a function and get the return value
get_free_card
echo  "cards present: $cards_present"


#set local variables
#EXE="mpirun -np 4 /home/clusterusers/jberryman/amber/bin/pmemd.MPI"
EXE=/home/clusterusers/jberryman/amber/bin/pmemd.cuda


sh ~/amber/amber.sh

from=crd
dest=heat
if [ -s "$stem.$dest" ]
then 
	echo "found $stem.$dest"
else
$EXE -O -i heat.in \
   -o $stem$dest.out -r $stem.$dest -x $stem$dest.nc \
   -c $stem.$from -ref $stem.crd -p $stem.top
fi

from=$dest
dest=heat2
if [ -s "$stem.$dest" ]
then
        echo "found $stem.$dest"
else
$EXE -O -i heat2.in \
   -o $stem$dest.out -r $stem.$dest -x $stem$dest.nc \
   -c $stem.$from -ref $stem.crd -p $stem.top
fi

from=$dest
dest=heat3
if [ -s "$stem.$dest" ]
then 
        echo "found $stem.$dest"
else
$EXE -O -i heat3.in \
   -o $stem$dest.out -r $stem.$dest -x $stem$dest.nc \
   -c $stem.$from -ref $stem.crd -p $stem.top
fi 

from=$dest
dest=heat4
if [ -s "$stem.$dest" ]
then
          echo "found $stem.$dest"
else
	  $EXE -O -i heat3.in\
	       -o $stem$dest.out -r $stem.$dest -x $stem$dest.nc \
	       -c $stem.$from -ref $stem.crd -p $stem.top
fi

from=$dest
dest=heat5
if [ -s "$stem.$dest" ]
then
          echo "found $stem.$dest"
else
          $EXE -O -i heat3.in\
               -o $stem$dest.out -r $stem.$dest -x $stem$dest.nc \
               -c $stem.$from -ref $stem.crd -p $stem.top
fi

from=$dest
dest=heat6
if [ -s "$stem.$dest" ]
then
          echo "found $stem.$dest"
else
          $EXE -O -i heat4.in\
               -o $stem$dest.out -r $stem.$dest -x $stem$dest.nc \
               -c $stem.$from -ref $stem.crd -p $stem.top
fi





zMax=$(tail -n 1 $stem.crd | awk '{printf("%s\n",substr($0,25,12))}')
for i in `seq 1 16`
do
	from=$dest
	dest=npt$i
	if [ -s "$stem.$dest" ]
	then
        	echo "found $stem.$dest"
	else

		echo "running $stem.$dest from $from"

$EXE -O -i npt.in \
   -o $stem$dest.out -r $stem.$dest -x $stem$dest.nc \
   -c $stem.$from -ref $stem.crd -p $stem.top

        ##reset the box length 
       boxLine=$(tail -n 1 $stem.$dest |\
       awk -v z="$zMax" '{printf("%s%s%s\n",substr($0,0,24),z,substr($0,37))}')
       N=$(wc -l $stem.$dest | awk '{print $1-1}')
       head -n $N $stem.$dest > c
       echo "$boxLine"      >> c
       echo "reset box from:"
       tail -n 1 $stem.$dest
       echo "to:"
       tail -n 1 c
       mv c $stem.$dest
        fi
done




running_pids=""

strtSnap=$dest
block_wait=0
for i in `seq 1 16`
do

        ##blocking wait for a non-busy card
	while [ "$block_wait" -eq "1" ] 
	do
		sleep 30
		count=0
		for pid in $running_pids
		do
		    c=$(ps -U jberryman | grep -c $pid)
		    count=$[c + $count]
		done
		if [ "$count" -lt "$count_cards_present" ]
		then
			##set active card and stop waiting
			get_free_card
			break
		fi
	done

##update the start snap to end of last equil... if not 
##running in parallel.
if [ -s "$stem.stretch${i}_0" ]
then
  strtSnap="stretch${i}_0"
fi

dest=$strtSnap
( for step in `seq 0 15`
do

from=$dest
dest=stretch${i}_$step
if [ -s "$stem.$dest" ] 
then
	echo "found $stem.$dest , skipping."
else

sed 's/t_stretch/t_stretch'$i'/'  stretch_$step.in > card$i.in

$EXE -O -i card$i.in \
   -o $stem$dest.out -r $stem.$dest -x $stem$dest.nc \
   -c $stem.$from -ref $stem.crd -p $stem.top 
  

fi
done ) &
pid=$!
block_wait=1
running_pids="$pid $running_pids"

done





