#! /bin/bash

progressdir="$1/progress"
logfile_name=$2
logfile="${progressdir}/${logfile_name}"

mkdir -p "$progressdir"

init_logging () {
    processname=$1
    touch "$logfile"
    #count=$(grep -c "$processname" "$logfile")
    #(( count+=1 ))
    {
        printf "\n%*s\n" 80 ' ' | tr ' ' '%'
        # printf "STARTING %s %s - %s" "$processname" "$count" "$(date '+%Y-%m-%d %H:%M:%S')"
        printf "STARTING %s - %s" "$processname" "$(date '+%Y-%m-%d %H:%M:%S')"
        printf "\n%*s\n\n" 80 ' ' | tr ' ' '%'
    } >> "$logfile"
}

finalize_logging () {
    processname=$1
    count=$(grep -c "$processname" "$logfile")
    (( count+=1 ))
    {
        printf "\n%*s\n" 80 ' ' | tr ' ' '-'
        # printf "%% FINISHED %s %s - %s" "$processname" "$count" "$(date '+%Y-%m-%d %H:%M:%S')"
        printf "%% FINISHED %s - %s" "$processname" "$(date '+%Y-%m-%d %H:%M:%S')"
        printf "\n%*s\n\n" 80 ' ' | tr ' ' '-'
    } >> "$logfile"
}

write_header () {
    header=$1
    prefix="$header "
    pad=$((80 - ${#prefix}))
    {
        printf "%s" "$prefix"
        printf "%*s\n\n" "$pad" '' | tr ' ' '-'
    } >> "$logfile"
}

write_progress () {
    printf "%s : %b\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$1" >> "$logfile"
}

skip () {
    write_progress "$1 already done -> skipping this step\n"
}

already_done () {
    [ -f "${progressdir}/$1.done" ]
}

eval_outcome () {
    if [ "$1" != 0 ]; then
        write_progress "Failed $2 -> exit value $1"
        exit 1
    else
        touch "${progressdir}/$2.done"
        write_progress "Finished $2 successfully\n"
    fi
}
