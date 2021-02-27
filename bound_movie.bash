#!/bin/bash

sys=$(pwd | awk -F '/' '{print $NF}')

for run in 1
do

rm -f bound_movie_$run.pdb

    for frame in `seq 1 3000`
    do

rm -f snapm.pdb
##unpack a pdb snapshot
(cpptraj ${sys}_noWat.pdb <<EOF
 trajin   ${sys}stretch_${run}_1to_15_noWat.nc $frame $frame 1
 strip    @BR*
 strip    @H*
 strip    @Na*
 strip    @F*
 trajout  snapm.pdb
EOF
) > ptraj.log

if [ -s "snapm.pdb" ]
then
echo $frame
else
break
fi

##filter out a pdb of docked ethidiums only.
##criterion is ET com distance from com of bps above and below
##also N23 and N24 distances from either COM of left and right base steps.
##
##hack to make a movie is print each residue twise, with resit ET and UT.
##ET == bonded (or all coords zeroed)
##UT == either bonded or not bonded
awk -v cut2=64 -v cutLoose=100 -v cutC2=64 '$1=="ATOM"{x=substr($0,30,9);\
                y=substr($0,39,8);\
                z=substr($0,47,9);i=int(substr($0,22,5));\
                line[i]=line[i]"\n"$0;
                bound[i]=0;
                comx[i]+=x;comy[i]+=y;comz[i]+=z;com_c[i]+=1}\
     $3=="N24"{x24[i]=x;y24[i]=y;z24[i]=z}\
     $3=="N23"{x23[i]=x;y23[i]=y;z23[i]=z}\
     $3=="O5'\''"{xp[i]=x;yp[i]=y;zp[i]=z}\
         END{for(p in xp){
               if(p%24==1||p%24==0){continue}\
               for(n in x24){\
                  dr2=(x24[n]-xp[p])**2+(y24[n]-yp[p])**2+(z24[n]-zp[p])**2;\
                  if(dr2>cutLoose){continue}\
                    pp=50-p;\
        dr2=(x23[n]-xp[pp])**2+(y23[n]-yp[pp])**2+(z23[n]-zp[pp])**2;\
                    if(dr2>cutLoose){continue}\
                    if(p>pp){q=pp;pp=50-q}else{q=p};\
                    comx_up=(comx[q]+comx[pp-1])/(2*com_c[q]);\
                    comx_dn=(comx[q-1]+comx[pp])/(2*com_c[q-1]);\
                    comy_up=(comy[q]+comy[pp-1])/(2*com_c[q]);\
                    comy_dn=(comy[q-1]+comy[pp])/(2*com_c[q-1]);\
                    comz_up=(comz[q]+comz[pp-1])/(2*com_c[q]);\
                    comz_dn=(comz[q-1]+comz[pp])/(2*com_c[q-1]);\
                    dc2 =((comx_up+comx_dn)*0.5-(comx[n]/com_c[n]))**2;\
                    dc2+=((comy_up+comy_dn)*0.5-(comy[n]/com_c[n]))**2;\
                    dc2+=((comz_up+comz_dn)*0.5-(comz[n]/com_c[n]))**2;\
                    if(dc2>cutC2){continue};\
                    comx_l=(comx[q]+comx[q-1])/(2*com_c[q]);\
                    comy_l=(comy[q]+comy[q-1])/(2*com_c[q]);\
                    comz_l=(comz[q]+comz[q-1])/(2*com_c[q]);\
                    dc24  =(comx_l-x24[n])**2;\
                    dc24 +=(comy_l-y24[n])**2;\
                    dc24 +=(comz_l-z24[n])**2;\
                    dc23  =(comx_l-x23[n])**2;\
                    dc23 +=(comy_l-y23[n])**2;\
                    dc23 +=(comz_l-z23[n])**2;\
                    if(dc24 > cut2 && dc23 > cut2){ continue }\
                    if(dc24<dc23){nn1=dc24}else{nn1=dc23}\
                    comx_r=(comx[pp-1]+comx[pp])/(2*com_c[pp-1]);\
                    comy_r=(comy[pp-1]+comy[pp])/(2*com_c[pp-1]);\
                    comz_r=(comz[pp-1]+comz[pp])/(2*com_c[pp-1]);\
                    dc24  =(comx_r-x24[n])**2;\
                    dc24 +=(comy_r-y24[n])**2;\
                    dc24 +=(comz_r-z24[n])**2;\
                    dc23  =(comx_r-x23[n])**2;\
                    dc23 +=(comy_r-y23[n])**2;\
                    dc23 +=(comz_r-z23[n])**2;\
                    if(dc24 > cut2 && dc23 > cut2){ continue }\
                    if(dc24<dc23){nn2=dc24}else{nn2=dc23}\
                    bound[n]=1;\
                    print n" is bound, with l,r,c,L: ",sqrt(nn1),sqrt(nn2),sqrt(dc2),sqrt(dr2);\
                    break;\
               }\
               print line[p];\
             }\
             for(n in x24){\
                print("checking resid ",n);\
                if(bound[n] == 0){\
                  print("not bound");\
                  xpl=1;out="";\
                  for(xpos=30;xpos<length(line[n]);xpos+=81){\
                    out = out""substr(line[n],xpl,30)"    0.000   0.000   0.000";xpl+=81;       }print out;\
                }else{print line[n]}\
                gsub(/ET/, "UT",line[n]);print line[n];\
             }\
            }' snapm.pdb   >> bound_movie_$run.pdb
            echo "ENDMDL" >> bound_movie_$run.pdb
    done
done
