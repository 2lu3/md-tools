#!/bin/bash

{% raw %}
pushd $(dirname $0)

jobs=($(ls job*.sh | sort -n))

# TODO: start, endなしに書き換える


if [ "$1" == "--local" ]; then
    for job in ${jobs[@]}; do
        echo "Running ${job}..."
        ./${job} --local
    done
else
    # if the len(jobs) == 1
    if [ ${#jobs[@]} -eq 1 ]; then
        pjsub ${jobs[0]}
    else
        # submit the first job
        msg=$(pjsub --step ${jobs[0]})
        jobid=$(echo $msg | cut -d " " -f6 | cut -d "_" -f1)
        # if there is at least one job
        if [ ${#jobs[@]} -gt 1 ]; then
            # submit the rest of the jobs
            for job in ${jobs[@]:1}; do
                echo "Submitting ${job}..."
                pjsub --step --sparam jid=$jobid ${job}
            done
        fi
    fi
fi
popd
{% endraw %}
