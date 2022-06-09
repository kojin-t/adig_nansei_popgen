#!/usr/bin/perl

open(LIST, $ARGV[0]);
while($l=<LIST>){
        chomp($l);
        @c=split(/\t/,$l);
        $a{$c[2]}=$l;
}
close LIST,

open(L, $ARGV[1]);
while($line=<L>){
        chomp($line);
        @b=split(/\t/,$line);
        if(exists($a{$b[2]})){
                print "$line\n";
        }
}
close(L);

