#!/usr/bin/awk -f


BEGIN { OFS="\t"; } {
  if ( $1 ~ />Cluster/ ) {
    ##################
    # print output
    ##################
    if ( NR > 1 ) {
      for ( s in abund_sample_array ) {
        abund = abund_sample_array[s]
        sample = sub("-Cluster.*","",s)
        print cluster,s,abund_sample_array[s],rep

      }
    }

    ####################
    # reset variables
    ####################
    abund=0
    delete abund_sample_array

    cluster=$0;
    sub(">","",cluster)

  } else {

    ####################
    # get representative
    ####################
     if ( $0 ~ "\*$" ) {
      rep=gensub(/\.\.\.$/,"","g",$3)
      rep=gensub(/>/,"","g",rep)
    }

    ####################
    # get abund per sample
    ####################
    sample=$3
    sub("_id.*","",sample)
    sub(">","",sample)
    abund_sample_array[sample]++

  }
}
END {

  for ( s in abund_sample_array ) {
        abund = abund_sample_array[s]
        sample = sub("-Cluster.*","",s)
        print cluster,s,abund_sample_array[s],rep
  }

}
