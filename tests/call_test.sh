../kled -R call_test.fa.gz -t 4 --NOF call_test.cram > call_test.vcf
diff <(sed '2d' call_test.answer.vcf) <(sed '2d' call_test.vcf)
