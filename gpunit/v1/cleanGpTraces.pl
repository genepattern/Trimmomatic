use strict;
use warnings;

while (my $line = <>) {
    next if ($line =~ /Process exited with status code:/);  # Local server executor appends this as final line

    $line =~ s/\[[\d]{2}:[\d]{2}:[\d]{2}\]/\[timestamp\]/g;
    $line =~ s/(\/{0,1})Axis[0-9]+\.att\_/$1/g;
    $line =~ s/\/[A-Za-z\/0-9@._-]+\//\[system_path\]\//g;
    
    print $line;
}
