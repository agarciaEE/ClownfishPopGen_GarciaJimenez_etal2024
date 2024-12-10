#!/bin/bash

n=$1  # Replace with your desired value for n
NR_IND=$2  # Replace with your desired value for NR_IND

INDS=()
for ((i = 0; i <= n; i += 2)); do
    INDS+=("$i")
done

SEL_INDS=()
for i in $(seq 1 $NR_IND); do
    rand=$((RANDOM % ${#INDS[@]}))
    SEL_INDS+=(${INDS[$rand]},$((INDS[$rand] + 1)))

    if [ $rand -eq 0 ]; then
        NEW_INDS=("${INDS[@]:(($rand+1))}")
    else
        NEW_INDS=("${INDS[@]:0:$((rand-1))}" "${INDS[@]:$rand}")
    fi

    INDS=()
    for e in "${NEW_INDS[@]}"; do
        if [ -n "$e" ]; then
            INDS+=("$e")
        fi
    done
done

INDEX=$(IFS=,; echo "${SEL_INDS[*]}")

echo "$INDEX"

