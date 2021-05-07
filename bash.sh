#passing arrays between functions
aa=({1..5})
bb=({A..E})

takeval() { # accept array
  i=$1[@]
  i=(${!i})
  j=$2[@]
  j=(${!j})
  echo ${i[*]} \| ${j[*]}
}

x=$(takeval aa bb) # convert back
a2=(${x%|*})     # cut out 1st array
b2=(${x#*|})
#echo ${aa[*]} ${a2[*]}
#echo ${bb[*]} ${b2[*]}
