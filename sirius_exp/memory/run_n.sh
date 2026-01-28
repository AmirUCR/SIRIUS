#!/bin/bash

mCitrine="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLMCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSKLSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
TIME=4800

counter=1
for i in 2 3 4 5 6 7 8 9 10; do
    outdir="sirius_out_${counter}"

    start=$(date +%s)

    # Run and capture peak RSS (kB) with GNU time
    peak_file=".peak_kb"
    /usr/bin/time -f "%M" -o "$peak_file" \
        sirius --prot "$mCitrine" -n $i --max_sec "$TIME"

    end=$(date +%s)

    dur=$((end - start))
    printf -v hms "%02d:%02d:%02d" $((dur/3600)) $(((dur%3600)/60)) $((dur%60))

    # Read peak memory and convert to MB
    peak_kb=$(tr -d '\n\r ' < "$peak_file")
    if [[ -z "$peak_kb" ]]; then peak_kb=0; fi
    peak_mb=$(( (peak_kb + 1023) / 1024 ))

    {
        echo "i=$i"
        echo "max_sec=$TIME"
        echo "start_epoch=$start"
        echo "end_epoch=$end"
        echo "duration_seconds=$dur"
        echo "duration_hms=$hms"
        echo "peak_memory_kb=$peak_kb"
        echo "peak_memory_mb=$peak_mb"
    } > "$outdir/measurements.txt"

    rm $peak_file

    ((counter++))
done