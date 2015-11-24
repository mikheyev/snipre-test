#!/bin/sh
# properties = {"output": ["../out/bayesian_results.csv"], "params": {}, "resources": {}, "input": ["../data/var/silentReplacement.csv", "../data/var/snps.csv", "../data/var/annotation.csv"], "rule": "snipre", "threads": 1, "local": false}
cd /work/MikheyevU/sasha/snipre_test/src && /home/s/sasha/sasha_env/bin/snakemake --snakefile /work/MikheyevU/sasha/snipre_test/src/Snakefile --force -j --keep-target-files --wait-for-files ../data/var/silentReplacement.csv ../data/var/snps.csv ../data/var/annotation.csv --latency-wait 5 --benchmark-repeats 1   --nocolor --notemp --quiet --no-hooks --nolock ../out/bayesian_results.csv --printshellcmds  --allowed-rules snipre  && touch "/work/MikheyevU/sasha/snipre_test/src/.snakemake/tmp.5PX3B9/4.jobfinished" || touch "/work/MikheyevU/sasha/snipre_test/src/.snakemake/tmp.5PX3B9/4.jobfailed"

