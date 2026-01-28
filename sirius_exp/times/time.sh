#!/bin/bash

SIRIUS=sirius

mCitrine="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLMCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSKLSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"

counter=1
for i in 5 5 5 10 10 10 20 20 20 40 40 40 80 80 80 160 160 160 320 320 320 640 640 640; do
    TIME=$((60 * i)) 
    outdir="sirius_out_${counter}"

    start=$(date +%s)
    $SIRIUS -prot $mCitrine -n 10 -t $TIME
    end=$(date +%s)

    dur=$((end - start))
    printf -v hms "%02d:%02d:%02d" $((dur/3600)) $(((dur%3600)/60)) $((dur%60))

    {
        echo "i=$i"
        echo "max_sec=$TIME"
        echo "start_epoch=$start"
        echo "end_epoch=$end"
        echo "duration_seconds=$dur"
        echo "duration_hms=$hms"
    } > "$outdir/runtime.txt"

    ((counter++))
done
