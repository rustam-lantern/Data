#!/bin/bash

FN=${1}
#declare -a s1=(15.498519 0.481036  0.000016097)
declare -a s1=(-3.511057 0.714451  0.000025758)
declare -a s2=(-2.463011 0.504433 0.0000079580)

declare -a s1=(15.498519 0.481036  0.000016097)
declare -a s2=(2.734449 0.485724 0.0000147857)


for (( i = 0; i < ${#s1[@]}; i++ ))
do
    if [[ $(more ${FN}) =~ ${s1[i]} ]] ; then
       #echo $(grep ${s1[i]} ${FN})
       sed -i "s/${s1[i]}/${s2[i]}/g" ${FN}
    fi
done

