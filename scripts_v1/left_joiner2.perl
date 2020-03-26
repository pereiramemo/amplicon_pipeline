#!/usr/bin/perl

# returns all rows from the left table (table1), with the matching rows \
# in the right table (table2). The result is NA in the right side when \
# there is no match.


open (FILE_1,$ARGV[0]) or die "Error: missing left table
usage:
perl left_joiner.perl table1 table2 column1 column2\n";

open(FILE_2,$ARGV[1]), or die "Error: missing right table
usage:
perl left_joiner.perl table1 table2 column1 column2\n";

$ARGV[2] or die "Error: missing right column number
usage:
perl left_joiner.perl table1 table2 column1 column2\n";

$ARGV[3] or die "Error: missing left column number
usage:
perl left_joiner.perl table1 table2 column1 column2\n";


$COLUMN_TABLE1=$ARGV[2];
$COLUMN_TABLE2=$ARGV[3];

while (<FILE_2>) {

    chomp $_;
    @line = split(/\t/,$_);
    $key = $line[ $COLUMN_TABLE2 - 1 ];
    push ( @{ $hash_table{ $key } }, $_ );
}



while (<FILE_1>) {
    
    
    chomp $_;
    @line = split(/\t/,$_);
    $key = $line[ $COLUMN_TABLE1 - 1 ];
    
    if ( @{ $hash_table{ $key } } ) {
     
      foreach $i ( @{$hash_table{ $key } } ) {
        
        @line2 = split(/\t/,$i);
        #$N = @line2 - 1;
        #delete $line2[ $COLUMN_TABLE2 -1 ];
        #splice(@line2,$COLUMN_TABLE2 -1,1);
        $line2 = join ("\t",@line2);
        print $_,"\t",$line2,"\n";
        
      }
      
    } else {
    
     # print $_,"\t","NA\t" x $N,"NA","\n";
     print $_,"\t","NA","\n";
    }  
    
}
