#!/bin/bash

sys=$(pwd | awk -F '/' '{print $NF}')

baseAtoms="N1,N2,N3,N4,N5,N6,N7,N8,N9"
baseAtoms="$baseAtoms,C1,C2,C3,C4,C5,C6,C7,C8,C9"
baseAtoms="$baseAtoms,O6"

keepThreshold=18.0

rm -f common_ats_2et.nc
rm -f common_ats_2et_furthest.nc
rm -f common_ats_2et_closest.nc
for run in `seq 1 16`
do

for bpStep in `seq 1 23`
do


	bp2=$[48 - $bpStep]

cat <<EOF > strip.pts
 parm ${sys}_noWat.pdb bondsearch 2.
 solvent  :ET
 trajin   ${sys}stretch_${run}_1to_15_noWat.nc 
 strip    @H*
 strip    @BR
 strip    @Na*
 strip    @Cl*
 strip    @F*
EOF


      if [ "$bpStep" -gt "1" ]
      then
	strip2_from=$[bp2 + 2]
cat <<EOF >> strip.pts
 strip    :${strip2_from}-48
EOF
      fi
      if [ "$bpStep" -lt "23" ]
      then
	strip2_to=$[bp2 - 1]
cat <<EOF >> strip.pts
 strip    :25-${strip2_to}
EOF
      fi
      if [ "$bpStep" -lt "23" ]
      then
	strip1_from=$[bpStep + 2]
cat <<EOF >> strip.pts
 strip    :${strip1_from}-24
EOF
      fi
      if [ "$bpStep" -gt "1" ]
      then
	strip1_to=$[bpStep - 1]
cat <<EOF >> strip.pts
 strip    :1-${strip1_to}
EOF
      fi

##unpack a traj: temp only, overwrite mode.
cat <<EOF >> strip.pts
 strip    :DG,DC,DA,DT,DG5,DG3,DC5,DC3@N9,C8,C7,N7,N6,O6,N4,O4,O2,N2,P,OP1,OP2
 center   :1-4
 trajout  step${bpStep}.nc
 trajout  step${bpStep}.pdb onlyframes 1
EOF
cpptraj  < strip.pts


cpptraj <<EOF
parm step${bpStep}.pdb bondsearch 2.
solvent  :ET
trajin   step${bpStep}.nc
closest  2 :1-4@$baseAtoms center outprefix onestep_2et closestout onestep_2et_stripping
trajout  common_ats_2et_closest.nc  box
trajout  common_ats_2et.pdb onlyframes 1
EOF

##strip everything where the ET is not inside keepThreshold from the 
##DNA
keepFrames=$(grep -v '#' onestep_2et_stripping |\
            awk '$3<='$keepThreshold'{keep[$1]=1}NF==4{max=$1}\
                                END{for(i=1;i<=max;i++){\
                                  if(i in keep){printf("%i ",i)}}\
                                   printf("\n")}')
stripFrames=$(grep -v '#' onestep_2et_stripping |\
            awk '$3<='$keepThreshold'{keep[$1]=1;total+=1}NF==4{max=$1;total+=1}\
                                END{for(i=1;i<=max;i++){\
                                  if(i in keep){n=0}else{printf("%i ",i)}}\
                                   printf("\n")}')

cpptraj <<EOF
parm    common_ats_2et.pdb bondsearch 2.
trajin  common_ats_2et_closest.nc onlyframes $keepFrames
trajout common_ats_2et.nc append
go
EOF


##for debug, keep a traj of frames that have been discarded
cpptraj <<EOF
parm    common_ats_2et.pdb bondsearch 2.
trajin  common_ats_2et_closest.nc onlyframes $stripFrames
trajout common_ats_2et_furthest.nc  append
go
EOF


rm  common_ats_2et_closest.nc
rm  step${bpStep}.nc
rm  step${bpStep}.pdb

    done

##tidy up.
rm -f step${bpStep}_1et.nc
rm -f bp.nc
done

##prepare some ad-hoc reaction coordinates
cpptraj <<EOF
parm    common_ats_2et.pdb bondsearch 2.
trajin  common_ats_2et.nc
distance :1-4@$baseAtoms :5 out dist_nuc_5.dat
distance :1-4@$baseAtoms :6 out dist_nuc_6.dat
distance :1,4@$baseAtoms :2,3@$baseAtoms out dist_rise.dat
go
EOF

cpptraj <<EOF
parm    common_ats_2et.pdb bondsearch 2.
trajin  common_ats_2et_furthest.nc
distance :1-4@$baseAtoms :5 out dist_nuc_52.dat
distance :1-4@$baseAtoms :6 out dist_nuc_62.dat
distance :1,4@$baseAtoms :2,3@$baseAtoms out dist_rise2.dat
go
EOF

cpptraj <<EOF
parm    common_ats_2et.pdb bondsearch 2.
trajin  common_ats_2et.nc
trajin  common_ats_2et_furthest.nc
rms first :1-4
principal :5 out princip_et1.dat
principal :6 out princip_et2.dat
vector    :1,3@$baseAtoms :2,4@$baseAtoms out dna_ax.dat
EOF


grep "ECTOR 2" princip_et1.dat | awk '{print $4,$5,$6}' > et1vec.dat
grep "ECTOR 2" princip_et2.dat | awk '{print $4,$5,$6}' > et2vec.dat
grep -v '#' dna_ax.dat         | awk '{r=sqrt($2*$2+$3*$3+$4*$4);print $2/r,$3/r,$4/r}' > dnavec.dat


paste dist_nuc_5.dat  dist_nuc_6.dat  | grep -v '#'  > dist56.dat
paste dist_nuc_52.dat dist_nuc_62.dat | grep -v '#' >> dist56.dat

grep -v '#' dist_rise.dat  >  c
grep -v '#' dist_rise2.dat >> c
rm   dist_rise2.dat
mv c dist_rise.dat


paste et1vec.dat et2vec.dat dnavec.dat dist56.dat dist_rise.dat |\
             awk '$11<$13{print $1*$7+$2*$8+$3*$9, $11, $15}\
                 $11>=$13{print $4*$7+$5*$8+$6*$9, $13, $15}' > angle_vs_distance.dat

cat angle_vs_distance.dat | awk '$3>5.&&$2<4.5{print NR,$1,$2,$3}' > possibles.dat

allowFrames=$(cat possibles.dat | awk '{printf("%i,",$1)}')
cpptraj <<EOF
parm    common_ats_2et.pdb bondsearch 2.
trajin  common_ats_2et.nc
trajin  common_ats_2et_furthest.nc
rms :1-4 first 
trajout possibles.nc onlyframes $allowFrames
EOF



##outputs a list of centre assignments mog_ts.dat and centroid frame id mog_labels.dat
../cluster_mog.py
label_max=$(awk 'BEGIN{m=-1}$1>m{m=$1}END{print m}' mog_labels.dat)
for i in `seq 0 $label_max`
do
   frames=$(awk '$1=='$i'{printf("%i,",NR)}' mog_labels.dat)
cpptraj <<EOF
parm    common_ats_2et.pdb bondsearch 2.
trajin  common_ats_2et.nc
trajin  common_ats_2et_furthest.nc
rms :1-4 first 
trajout centroid_$i.nc onlyframes $frames
EOF


done




exit

##do clustering
cpptraj common_ats.pdb <<EOF
trajin  common_ats.nc
parm common_ats.pdb  bondsearch 3.
center  :1-4
rms first
 cluster kmeans clusters 20 randompoint kseed 8465623 \
        sieve 40 random sieveseed 934478 \
        out cnumvtime.dat \
        summary summary.dat \
        info info.dat    \
        cpopvtime cpopvtime.agr normframe \
        repout rep repfmt pdb \
        singlerepout singlerep.nc singlerepfmt netcdf \
        avgout Avg avgfmt restart
EOF

##do clustering
cpptraj common_ats_2et.pdb <<EOF
trajin  common_ats_2et.nc
parm common_ats_2et.pdb  bondsearch 3.
center  :1-4
rms first
 cluster kmeans clusters 20 randompoint kseed 8465623 \
        sieve 40 random sieveseed 934478 \
        out cnumvtime_2et.dat \
        summary summary_2et.dat \
        info info_2et.dat    \
        repout rep_2et repfmt pdb \
        singlerepout singlerep_2et.nc singlerepfmt netcdf \
        avgout Avg_2et avgfmt restart
EOF

##inspect and set by hand: vmd rep.c[0-9].pdb rep.c[0-9][0-9].pdb

if [ "$sys" == "gggggg_et" ] 
then
    interced_clusters="0 3 5"
elif [ "$sys" == "ggcgac_et" ] 
    interced_clusters="0 4 8 15"
fi

rm clusCount_*.dat
    
   grep -v "\#" cnumvtime_2et.dat |\
        awk -v interc="$interced_clusters"\
         'BEGIN{split(interc, a)}\
               {t=(($1-1)%3000)+1;isin=0;clus=$2;\
                   for(c in a){\
                               if(a[c] == clus){isin=1;break}}
                   if(isin==1){count_i[t]+=1}\
                          else{count_o[t]+=1}\
                   count[clus][t]+=1}\
                  END{printf("# ");\
                       for(c in a){printf("%i ",a[c])};\
                       printf("\n");\
                       for(t=1;t<=3000;t++){\
                       if(t in count_i || t in count_o){\
                          printf("%i ",t);\
                          for(c in a){\
                            printf("%.4f ",count[a[c]][t]/(16.*23))\
                          }\
                          if(t in count_i){\
                printf("%.4f ", count_i[t]/(16.*23))}else{printf("0.0000 ")}\
                          if(t in count_o){\
                printf("%.4f\n",count_o[t]/(16.*23))}else{printf("0.0000\n")}\
                        }}}' > clusCount_all.dat
