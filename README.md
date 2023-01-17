# lrdf
Long Read (ONT) Data Frame

Install:

```
$ cargo install lrdf
```

Run:

```
$ lrdf example-ont.fastq < example-ont.fastq.lrdf.tsv
```

```
$ head 20221114_DNA_Wolf_Kemter_7plex-BC06.porechop.fastq.lrdf.tsv
seq_length mean_quality kmer_start kmer_end nt_A nt_G nt_T nt_C nt_U channels start_times
3871       32.94        ACTA       AGAA     1275 767  1010 819  0    2272     2022-11-14 13:43:24 +00:00
1972       33.18        CCTG       CCAA     472  395  599  506  0    1229     2022-11-14 13:43:24 +00:00
336        37.12        TTAT       GGTA     92   47   131  66   0    951      2022-11-14 13:43:25 +00:00
3632       30.80        TTAT       CCAC     991  849  970  822  0    1182     2022-11-14 13:43:25 +00:00
3789       23.81        GGCA       AAAG     877  916  795  1201 0    1273     2022-11-14 13:43:25 +00:00
3364       32.50        AGTT       AGGA     1457 591  626  690  0    2955     2022-11-14 13:43:25 +00:00
4336       28.98        GGAA       TTCA     1472 713  1343 808  0    1214     2022-11-14 13:43:25 +00:00
2483       19.93        TTTT       TAAC     681  543  691  568  0    1594     2022-11-14 13:43:26 +00:00
4035       28.10        TTCA       GCCC     558  1705 602  1170 0    307      2022-11-14 13:43:26 +00:00
```
