viral_genus=(10804.*10803
1089136.*1623295
11053.*11051
1136133.*194960
11801.*153135
1327964.*1925779
1544901.*11102
198503.*10509
37138.*10806
915425.*325455)

viral_total=(taxid_10804
taxid_1089136
taxid_11053
taxid_1136133
taxid_11801
taxid_1327964
taxid_1544901
taxid_198503
taxid_37138
taxid_915425)

bacteria_genus=(1036673.*44249
1042876.*286
1053692.*2184
1104326.*547
1117943.*28105
1123863.*53335
1133568.*29521
1159202.*2093
1184253.*768
1218356.*810
122586.*482
1328311.*64895
226185.*1350
243276.*157
262728.*724
290339.*413496
299768.*1301
300852.*270
347495.*1386
354242.*194
367737.*28196
373384.*620
395492.*379
402882.*22
418127.*1279
434271.*713
440085.*407
441952.*262
449216.*780
536232.*1485
552536.*1637
568706.*517
573059.*872
573236.*1678
768492.*613
862962.*816
889738.*469
930943.*2284
990315.*338
998820.*1578)

bacteria_total=(taxid_1036673
taxid_1042876
taxid_1053692
taxid_1104326
taxid_1117943
taxid_1123863
taxid_1133568
taxid_1159202
taxid_1184253
taxid_1218356
taxid_122586
taxid_1328311
taxid_226185
taxid_243276
taxid_262728
taxid_290339
taxid_299768
taxid_300852
taxid_347495
taxid_354242
taxid_367737
taxid_373384
taxid_395492
taxid_402882
taxid_418127
taxid_434271
taxid_440085
taxid_441952
taxid_449216
taxid_536232
taxid_552536
taxid_568706
taxid_573059
taxid_573236
taxid_768492
taxid_862962
taxid_889738
taxid_930943
taxid_990315
taxid_998820)

result_file=$1

b_sum=0
for ss in ${bacteria_genus[@]}
do 
echo $ss
total=`grep $ss $result_file | wc -l`
echo $total
b_sum=`expr $b_sum + $total`
done

b_total=0
for ss in ${bacteria_total[@]}
do 
echo $ss
total=`grep $ss $result_file | wc -l`
echo $total
b_total=`expr $b_total + $total`
done

v_sum=0
for ss in ${viral_genus[@]}
do 
echo $ss
total=`grep $ss $result_file | wc -l`
echo $total
v_sum=`expr $v_sum + $total`
done

v_total=0
for ss in ${viral_total[@]}
do 
echo $ss
total=`grep $ss $result_file | wc -l`
echo $total
v_total=`expr $v_total + $total`
done

echo "viral total:" $v_total
echo "bacteria total:" $b_total

echo "viral sum:" $v_sum
echo "bacteria sum:" $b_sum
