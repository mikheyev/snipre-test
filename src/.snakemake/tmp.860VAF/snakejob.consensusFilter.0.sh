#!/bin/sh
# properties = {"output": ["../data/var/harpur.filtered.vcf"], "params": {}, "resources": {}, "input": ["../data/var/harpur.vcf.gz"], "rule": "consensusFilter", "threads": 1, "local": false}
cd /work/MikheyevU/sasha/snipre_test/src && /home/s/sasha/sasha_env/bin/snakemake --snakefile /work/MikheyevU/sasha/snipre_test/src/Snakefile --force -j --keep-target-files --wait-for-files ../data/var/harpur.vcf.gz --latency-wait 5 --benchmark-repeats 1   --nocolor --notemp --quiet --no-hooks --nolock ../data/var/harpur.filtered.vcf --printshellcmds  --allowed-rules consensusFilter  && touch "/work/MikheyevU/sasha/snipre_test/src/.snakemake/tmp.860VAF/0.jobfinished" || touch "/work/MikheyevU/sasha/snipre_test/src/.snakemake/tmp.860VAF/0.jobfailed"

