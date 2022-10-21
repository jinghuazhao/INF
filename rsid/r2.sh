#!/usr/bin/bash

export prot=VEGF.A
export trait=ieu-a-996
export rt=${INF}/mr/gsmr
export gsmr=${rt}/gsmr-5e-8-3/out

# GSMR
cat <(echo snp A1 A2 eaf bzx bzx_se bzy bzy_se) \
    <(gunzip -c ${gsmr}/${prot}-${trait}.eff_plot.gz | awk '/effect_begin/,/effect_end/' | grep -v effect) | sort -k1,1

# Protein
gunzip -c ${rt}/prot/${prot}.gz | \
grep -f <(gunzip -c ${gsmr}/${prot}-${trait}.eff_plot.gz | awk '/effect_begin/,/effect_end/' | grep -v effect | cut -d' ' -f1) | \
sort -k1,1
grep -f <(gunzip -c ${gsmr}/${prot}-${trait}.eff_plot.gz | awk '/effect_begin/,/effect_end/' | grep -v effect | cut -d' ' -f1) \
        ${rt}/prot/${prot}-rsid.txt | sort -k1,1

# Trait
gunzip -c ${rt}/trait/${prot}-${trait}.gz | \
grep -f <(gunzip -c ${gsmr}/${prot}-${trait}.eff_plot.gz | awk '/effect_begin/,/effect_end/' | grep -v effect | cut -d' ' -f1) | \
sort -k1,1
