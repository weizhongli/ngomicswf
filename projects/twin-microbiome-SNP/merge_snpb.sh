mkdir merge_snpb

for tid in `ls -1 */SNPb1/[1-9]*.indel2.tsv |  gawk -F '/' '{print $3}' | sed 's/.indel2.tsv//' | sort | uniq`; do
  # ls -1 */SNPb/$tid.vcf.tsv
  echo "/home/wli/git/ngomicswf/NGS-tools/NGS-table-merge.pl -s NGS-samples      -t 0 -p 0 -f SNPb1/$tid.indel2.tsv -k 1 -i 0 -v 1 -o merge_snpb/$tid.indel2.tsv"
  echo "/home/wli/git/ngomicswf/NGS-tools/NGS-table-merge.pl -s NGS-samples      -t 0 -p 0 -f SNPb1/$tid.vcf2.tsv   -k 1 -i 0 -v 1 -o merge_snpb/$tid.vcf2.tsv"

  /home/wli/git/ngomicswf/NGS-tools/NGS-table-merge.pl -s NGS-samples      -t 0 -p 0 -f SNPb1/$tid.indel2.tsv     -k 1 -i 0 -v 1 -o merge_snpb/$tid.indel2.tsv &
  /home/wli/git/ngomicswf/NGS-tools/NGS-table-merge.pl -s NGS-samples      -t 0 -p 0 -f SNPb1/$tid.vcf2.tsv       -k 1 -i 0 -v 1 -o merge_snpb/$tid.vcf2.tsv &
  wait
done


