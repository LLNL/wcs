#!/bin/sh

#seeds="59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 233 239 241 251 257 263 269 271 277 281 283 293 307 311 313 317 331 337 347 349 353 359 367 373 379 383 389 397 401 409 419 421 431 433 439 443 449 457 461 463 467 479 487 491 499 503 509 521 523 541 547 557 563 569 571 577 587 593 599 601 607 613 617 619 631 641 643 647 653 659 661 673 677 683 691 701 709 719 727 733 739 743 751 757 761 769 773 787 797 809 811 821 823 827 829 839 853 857 859 863 877 881 883 887 907 911 919 929 937 941 947 953 967 971 977 983 991 997"
#seeds="59 67 73 83 97 103 109 127 137 149 157 167 179 191 197 211 227 233 241 257"
seeds="59 67 73 83 97 103 109 127 137 149 157"


n=`echo $seed | wc -w`

#prob=eq38
prob=eq29

if [ "${prob}" == "eq29" ] ; then
  tmax=5
  sample_interval=0.1
  Y0="10 1000 3000"
elif [ "${prob}" == "eq38" ] ; then
  tmax=30
  sample_interval=0.5
  Y0="1000"
else
  echo "Unknown problem"
  exit
fi


odir=result-sample-${sample_interval}
mkdir -p ${odir}/toplot


for Y in ${Y0}
do
  for m in 0 1
  do
    # Write the header into the summary file,
    # which is the sequence of sampling times
    echo ${sample_interval} ${tmax} \
      | awk '{
          printf("0");
          for (t=$1; t <= $2; t+=$1) {
            printf(" %f", t);
          }
          printf("\n"); }' \
      > ${odir}/toplot/summary.${prob}.Y-${Y}.m-${m}.txt

    for s in ${seeds}
    do
      echo Y0 ${Y} method ${m} seed ${s}
      input=${prob}.${Y}.graphml
      output=${prob}.Y-${Y}.m-${m}.seed-${s}.txt
      cat ${prob}.graphml | sed -e 's/\@Y0\@/'$Y'/g' > ${odir}/${input}
      cmd_part1="bin/ssa -r t${sample_interval} -t $tmax -s $s -m $m "
      cmd_part2="${odir}/${input} -o ${odir}/${output}"
      cmd="${cmd_part1} ${cmd_part2}"
      echo ${cmd} > /dev/stderr
      ${cmd}

      if [ "$prob" == "eq29" ] ; then
        # t, Y
        awk '(FNR>2){print $1, $2}' ${odir}/${output} \
          > ${odir}/toplot/${output}
      elif [ "$prob" == "eq38" ] ; then
echo
        # t, Y1, Y2
        awk '(FNR>3){print $1, $2, $3}' ${odir}/${output} \
          > ${odir}/toplot/${output}
      fi

      # Summary of Y or Y1 population
      if [ -f ${odir}/toplot/${output} ] && [ ! -d ${odir}/toplot/${output} ]
      then
        data=`cat ${odir}/toplot/${output} | awk '{print $2}'`
        echo ${s} ${data} >> ${odir}/toplot/summary.${prob}.Y1-${Y}.m-${m}.txt
        data=`cat ${odir}/toplot/${output} | awk '{print $3}'`
        echo ${s} ${data} >> ${odir}/toplot/summary.${prob}.Y2-${Y}.m-${m}.txt
      fi
    done > ${odir}/results.${prob}.Y-${Y}.m-${m}.txt
  done
done
