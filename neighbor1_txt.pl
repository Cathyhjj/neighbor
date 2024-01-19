#! /usr/bin/perl -w

use strict;
use File::Glob ':globally';
use File::Basename;
use File::Copy;

my($i,$j,$k);
my($line);
my($name);
my(@x,@y,@z,@ipot,@tag,@dis,@atom); #original cluster group
my(@elmt);
my(@xs,@ys,@zs,@ipots);  #smaller cluster group
my($tagforipot1, $tagforipot2, $tagforipot3);
my(@r1,@r2,@r3,@n1,@n2,@n3);
my(%dis1,%dis2, %dis3);
my($center, $neighbor,$number_of_atoms);
my($total_ipot1);
my($a,$b,$c);
my(%x, %y,%z);
my(@x2d, @y2d, @z2d);
my(@x_planes, @y_planes, @z_planes);
my(%a_b);
my(@ipot1_atomnumbers);
my($key);

$name = $ARGV[0];
$atom[1] = $ARGV[1], if (defined($ARGV[1])); 
$atom[2] = $ARGV[2], if (defined($ARGV[2])); 
$atom[3] = $ARGV[3], if (defined($ARGV[3])); 

print "\nneighbor:  reading file $name\n";
open(DATA,$name) || die "cannot open $name";
while ($line = <DATA>){
    if ($line =~ /Elmt/){
	$line =<DATA>;  #skip next line of dashes
	$i=0;
	while ($line =<DATA>){
	    last, if($line =~/^-/);
	    $line =~ s/^[\s]+|[\s]+$//g;
	    $line =~ s/[\s]+/ /g;
	    #print "line: $line\n";
	    ($elmt[$i], $tag[$i], $x[$i],$y[$i],$z[$i]) = split(/\s/,$line);
	    #print "tag = $tag[$i], $i, ";
	    if ($tag[$i] eq $tag[0]){
		$tagforipot1 = $tag[0];
		$ipot[$i] = 1;
		#print "ipot1\n";
		$atom[1] = $elmt[0],  if (!(defined($atom['1']))); 	
	    }elsif(!defined($tagforipot2)){
		$tagforipot2 = $tag[$i];
		$ipot[$i] = 2;
		#print "ipot2\n";
		$atom[2] = $elmt[$i],  if (!(defined($atom['2']))); 
	    }elsif($tag[$i] eq $tagforipot2){
		$ipot[$i] = 2;
		#print "ipot2\n";
	    }elsif(!defined($tagforipot3)){
		$tagforipot3 = $tag[$i];
		$ipot[$i] = 3;
		#print "ipot3\n";
		$atom[3] = $elmt[$i],  if (!(defined($atom['3']))); 
	    }elsif($tag[$i] eq $tagforipot3){
		$ipot[$i] = 3;
		#print "ipot3\n";
	    }
	    $i++;
	}
    }
}

close(DATA);
$number_of_atoms=$i;
print "\tnumber of atoms in total cluster: $number_of_atoms\n";

my($repeat,@remove);
$repeat =1;  #do it again?
$remove[0] =0;  #remove some atoms from cluster

while ($repeat){
  print "repeat is true\n";
  #put original x,y,z,ipot into working groups
  my(@temparray);
  @temparray = &removeatoms2(\@remove,\%a_b, \@r1, \@x,\@y,\@z,\@ipot);
  #reset all variables 
  @r1 = 0;
  @r2 = 0;
  @r3 = 0;
  @n1 = 0;
  @n2 = 0;
  @n3 = 0;
  @xs = 0;
  @ys = 0;
  @zs = 0;
  @x2d = 0;
  @y2d = 0;
  @z2d = 0;
  @x_planes = 0;
  @y_planes = 0;
  @z_planes = 0;
  
  $total_ipot1=0;
  @ipot1_atomnumbers=0;
  @dis = 0;
  #number of atoms for each distance
  foreach $key (keys(%dis1)){
    delete ($dis1{$key});
  }
  foreach $key (keys(%dis2)){
    delete ($dis2{$key});
  }
  foreach $key (keys(%dis3)){
    delete ($dis3{$key});
  }
  #number of atoms for a given plane 
  foreach $key (keys %x){
    delete $x{$key};
  }
  foreach $key (keys %y){
    delete $y{$key};
  }
  foreach $key (keys %z){
    delete $z{$key};
  }

  #distance between a and b  key is make of atom numbers
  foreach $key (keys(%a_b)){
    delete ($a_b{$key});
  }
  @xs = &part_array(1,4,@temparray);
  @ys = &part_array(2,4,@temparray);
  @zs = &part_array(3,4,@temparray);
  @ipots = &part_array(4,4,@temparray);
  $number_of_atoms = @xs;
  print "\tnumber of atoms in smaller cluster is $number_of_atoms\n";

  #count number of atoms for each x, y, and z position
  ($i, $j, $k) = (0,0,0);
  my($atom,$key);
  for ($atom=0; $atom<$number_of_atoms; $atom++){
    $x2d[$atom] = &round($xs[$atom]);
    $y2d[$atom] = &round($ys[$atom]);
    $z2d[$atom] = &round($zs[$atom]);
    $x{$x2d[$atom]}++;
    $y{$y2d[$atom]}++;
    $z{$z2d[$atom]}++;
  }
  #put the keys into an array to be ordered
  foreach $key (sort(keys(%x))){
    $x_planes[$i++] = $key;
  }
  foreach $key (sort(keys(%y))){
    $y_planes[$j++] = $key;
  }
  foreach $key (sort(keys(%z))){
    $z_planes[$k++] = $key;
  }
  @x_planes[0..$i-1] = &order(@x_planes);
  @y_planes[0..$j-1] = &order(@y_planes);
  @z_planes[0..$k-1] = &order(@z_planes);
  
  #print the planes and the number of atoms on each plane
  print "x planes followed by number of atoms\n";
  for($i=0; $i < @x_planes; $i++){
    print "$x_planes[$i]\t";
  }
  print "\n";
  for($i=0; $i < @x_planes; $i++){
    print "$x{$x_planes[$i]}\t";
  }
  print "\n";

  print "y planes followed by number of atoms\n";
  for($i=0; $i < @y_planes; $i++){
    print "$y_planes[$i]\t";
  }
  print "\n";
  for($i=0; $i < @y_planes; $i++){
    print "$y{$y_planes[$i]}\t";
  }
  print "\n";
  print "z planes followed by number of atoms\n";
  for($i=0; $i < @z_planes; $i++){
    print "$z_planes[$i]\t";
  }
  print "\n";
  for($i=0; $i < @z_planes; $i++){
    print "$z{$z_planes[$i]}\t";
  }
  print "\n";

  #calculate distances between atoms
  $i=0;
  for ($center=0; $center<$number_of_atoms; $center++){
    next, unless (($ipots[$center] eq 1) || ($ipots[$center] eq 0));
    $ipot1_atomnumbers[$total_ipot1] = $center;
    $total_ipot1++;
    for($neighbor=0; $neighbor<$number_of_atoms; $neighbor++){
     #print "$neighbor, $center\n";
      next, if ($center =~ /^$neighbor$/);
      $dis[$i] = sqrt( ($xs[$center]-$xs[$neighbor])**2 +
		       ($ys[$center]-$ys[$neighbor])**2 +
		       ($zs[$center]-$zs[$neighbor])**2);
      $dis[$i] = &round($dis[$i]);
      $a_b{$center . "_" . $neighbor} = $dis[$i];
      $a_b{$center . "_" . $dis[$i]}++;
      $dis1{$dis[$i]}++, if(($ipots[$neighbor] eq 1) || ($ipots[$neighbor] eq 0));
      $dis2{$dis[$i]}++, if($ipots[$neighbor] eq 2);
      $dis3{$dis[$i]}++, if($ipots[$neighbor] eq 3);
      ##print "distance= $dis[$i],  $ipots[$neighbor] $neighbor\n";
      $i++;
    }
  }
  #calculate average number of neighbors and sort them;
  my($key1,$key2, $key3);
  my($k)=0;
  #print "calculating average\n";
  foreach $key1 (sort(keys(%dis1))){
     ($r1[$k], $n1[$k])=  ($key1, $dis1{$key1}/($total_ipot1));
     $k++;
   }
  my(@array,$length,$half);
  (@array) = &clean_array(@r1,@n1);
  (@r1,@n1) = (0,0);
  @r1 = &firsthalf(@array);
  @n1 = &secondhalf(@array);
  $k = @r1;

  my ($j)=0;
    foreach $key2 (sort(keys(%dis2))){
      ($r2[$j], $n2[$j]) = ($key2, $dis2{$key2}/($total_ipot1));
      #print "found distances\n";
      $j++;
    }
    (@array) = &clean_array(@r2,@n2);
    (@r2,@n2) = (0,0);
    @r2 =&firsthalf(@array);
    @n2 = &secondhalf(@array);
    $j = @r2;

  my ($l)=0;
    foreach $key3 (sort(keys(%dis3))){
      ($r3[$l], $n3[$l]) = ($key3, $dis3{$key3}/($total_ipot1));
      $l++;
    }
    (@array) = &clean_array(@r3,@n3);
    (@r3,@n3) = (0,0);
    @r3=&firsthalf(@array);
    @n3=&secondhalf(@array);
    $l = @r3;

  #print out a square block padded with zeros of the results
  my($max);
  $max=$k, if($j<$k && ($l<$k) );
  $max=$j, if($k<=$j && ($l<$j) );
  $max=$l, if($k<=$l && ($j<=$l) );
  print "#pot 1, 2  and 3\n#number and distances(Angstroms)\n";
  $atom['1'] = '', if (!(defined($atom['1'])));
  $atom['2'] = '', if (!(defined($atom['2'])));
  $atom['3'] = '', if (!(defined($atom['3'])));
  print "#N_$atom['1']$atom['1']\tR\tN_$atom['1']$atom['2']\tR\tN_$atom['1']$atom['3']\tR\n";
  for($i=0; $i<$max; $i++){
    ($r1[$i],$n1[$i])=(0,0), if(!(defined($n1[$i])));
    ($r2[$i],$n2[$i])=(0,0), if(!(defined($n2[$i])));
    ($r3[$i],$n3[$i])=(0,0), if(!(defined($n3[$i])));
    printf "%5.2f %5.2f   %5.2f %5.2f   %5.2f %5.2f\n",
      $n1[$i], $r1[$i], $n2[$i], $r2[$i], $n3[$i], $r3[$i];
  }


  #print out the number of atoms on surface acording to first neighbor coordination number
  print "near neighbor distance is $r1[0]\n";
  my(@neighbors)=0;
  my($temp);
  for ($i=0; $i<@ipot1_atomnumbers; $i++){
    $key = $ipot1_atomnumbers[$i] . "_" . $r1[0];
    if (defined($a_b{$key})){
      $temp = $a_b{$key};
      $neighbors[$temp]++;
    }
  }
  for ($i=0; $i<@neighbors; $i++){
    if(defined($neighbors[$i])){
      print "$neighbors[$i] atoms have $i neighbors\n";;
    }
  }
      
  @remove = 0;
  print "do you want to remove some atoms? Nx, Ny, or Nz\n";
  print "Nx, Ny, Nz to remove all negative x, y, or z atoms\n";
  print "CN # to remove all atoms with CN less than some #\n";
  print "R # to remove all atoms with radial distance greater than R from center\n";
  print "X|Y|Z # to remove all atoms on x|y|z plane with distance greater than # \n";
  $key = <STDIN>;
  if ($key =~ /^N/){
    print "found N removing negatitive values\n";
    $remove[0] = $key;
    $repeat = 1;
  }
  if ($key =~ /CN[\s]+(\d\d?)/){
    $remove[1] = $1;
    print "found CN\n";
    print "removing all surface atoms with CN less than $remove[1]\n";
    $repeat = 1;
  }
  if ($key =~ /R[\s]+(\d\d?.*)/){
    $remove[2] = $1;
    print "found Radial distance\n";
    print "removing all atoms with radial distance greater than  $remove[2]\n";
    $repeat = 1;
  }
  if ($key =~ /(X|Y|Z)[\s]+(\d\d?.*)/){
    $remove[3] = $1;
    $remove[4] = $2;
    print "found plane X|Y|Z  distance\n";
    print "removing all atoms on plane $remove[3] with radial distance greater than  $remove[4]\n";
    $repeat = 1;
  }
#  if ($key =~ /print/){
#    print "Enter file name:  ";
#    $name = <STDIN>;
#    open (OUTPUT,  "> $name") || die "Couldn't create output file\n";
#    my($num);
#    $num = @xs;
#    printf OUTPUT "%4d ATOMS,    0 BONDS,    0 CHARGES\n",$num;
#    for ($i=0; $i<@xs; $i++){
#      printf OUTPUT "%5d Re  %10.5f%10.5f%10.5f 0.00\n", $i, $xs[$i], $ys[$i], $zs[$i];
#    }
#    $repeat = 1;
#    close (OUTPUT);
#  }
  if ($key =~ /print/){
    print "Enter file name:  ";
    $name = <STDIN>;
    $name = $name . ".chem3d";
    #open (OUTPUT,  "> $name") || die "Couldn't create output file\n";
    my($num);
    $num = @xs;
    #printf OUTPUT "%4d\n",$num;
    printf "%4d\n",$num;
    for ($i=0; $i<@xs; $i++){
      if($ipots[$i] =~ /2/){
	printf "S %5d %10.5f%10.5f%10.5f\n", $i, $xs[$i], $ys[$i], $zs[$i];
      }
      if($ipot[$i] =~ /1/){
	printf "W %5d %10.5f%10.5f%10.5f\n", $i, $xs[$i], $ys[$i], $zs[$i];
      }
      if($ipot[$i] =~ /0/){
	printf "W %5d %10.5f%10.5f%10.5f\n", $i, $xs[$i], $ys[$i], $zs[$i];
      }
    }
    $repeat = 1;
    close (OUTPUT);
  }
  if ($key =~ /^$/){
    $repeat = 0;
  }
  print "value for repeat is $repeat\n";
}





###########################################################################
sub round{
  my($number) = @_;
  my $third_digit;
  #print "number: $number\n";
  $number =~ /(-?\d?\d?\.?\d?\d?)(\d?)/;
  $number = $1;
  $third_digit= $2;
  $third_digit= 0, if ($third_digit =~ //);
  #print "before $number\n";
  $number = $number + 0.01, if ($third_digit >= 5);
  #print "third $third_digit\n";
  #print "after $number\n";
  return($number);
}

sub clean_array{
  #sort will put numbers in order but 10's are first this sub fixes that
  my(@array)=@_;
  my($i,@r,@n);
  my(@rshort,@nshort);
  my($k,$j,$l,$element);
  $l = @array;
  $k = $l/2;
  @r =@array[0..$k-1];
  @n=@array[$k..$l];
  $r[0] = 0, if (!(defined($r[0])));
  $j=0;
  for($i=0; $i<$k; $i++){
    if($r[$i] < 10){
      $rshort[$j] = $r[$i];
      $nshort[$j] = $n[$i];
      $j++;
    }
  }
  return(@rshort,@nshort);
}

sub order{
  my(@array)=@_;
  my($i,$temp,$length);
  $length = @array;
  for ($i=0; $i<$length-1; $i++){
    $a = $array[$i] + 0;
    $b = $array[$i+1] + 0;
    if ($a < $b){
      $temp = $array[$i];
      $array[$i] = $array[$i+1];
      $array[$i+1] = $temp;
      $i=-1;
    }
  }
  return(@array);
}


sub part_array{
  my($section,$parts, @array)=@_;
  my($length, $newlength,@newarray);
  my($start,$end);
  $length = @array;
  $start = ((($section-1)*$length)/$parts);
  $start = 0, if ($section == 1);
  $end = ($section*$length)/$parts-1;
  $newlength = $end-$start+1;
  @newarray[0..$newlength-1] = @array[$start..$end];
  #print"\n\n";
  #print "section=$section  parts=$parts length=$length newlength=$newlength\n";
  #print "start=$start  end=$end\n";
  #print "@array";
  #print "\n\n";
  #print"@newarray";
  #print"\n\n";
  return(@newarray);
}


sub firsthalf{
  my(@array)=@_;
  my($length, $half,@firsthalf);
  $length = @array;
  $half = $length/2;
  @firsthalf[0..$half-1] = @array[0..$half-1];
  #print "first half\n";
  #print "@firsthalf";
  #print"\n";
  return(@firsthalf);
}

sub secondhalf{
  my(@array)=@_;
  my($length, $half,@secondhalf);
  $length = @array;
  $half = $length/2;
  @secondhalf[0..$half-1] = @array[$half..$length];
  #print "secondhalf\n";
  #print "@secondhalf";
  return(@secondhalf);
}




sub removeatoms2{
  my($i,$length);
  my($refx,$refy,$refz,$refipot);
  my($refa_b,$refkey,$r1);
  my($x,$y,$z,$ipot);
  my(@temp,$temp,$numberatoms);
  ($refkey, $refa_b, $r1, $refx, $refy, $refz, $refipot)=@_;
  my($j)=0;
  $length = @$refkey;
  $numberatoms = @$refx;
  print "there are $length things to do\n";
  if(defined($$refkey[1])){
    print "print removing CN atoms\n";
    for ($i=0; $i<$numberatoms; $i++){
      if ($$refa_b{$i . "_" . $$r1[0]} < $$refkey[1]){
	splice(@$refx, $j, 1);
	splice(@$refy, $j, 1);
	splice(@$refz, $j, 1);
	splice(@$refipot, $j, 1);
	$j=$j-1;
      }
      $j++;
    }
  }
  if(defined($$refkey[2])){
    print "print removing atoms with radial distance larger than $$refkey[2] \n";
    for ($i=0; $i<$numberatoms; $i++){
      if ($$refa_b{0 . "_" . $i} > $$refkey[2]){
	splice(@$refx, $j, 1);
	splice(@$refy, $j, 1);
	splice(@$refz, $j, 1);
	splice(@$refipot, $j, 1);
	$j=$j-1;
      }
      $j++;
    }
  }
  if(defined($$refkey[3])){
    print "print removing atoms on plane $$refkey[3] with radial distance larger than $$refkey[4] \n";
    if ($$refkey[3] =~ /X/){
      for ($i=0; $i<@$refx; $i++){
	if (abs($$refx[$i]) > $$refkey[4]){
	  splice(@$refx, $i, 1);
	  splice(@$refy, $i, 1);
	  splice(@$refz, $i, 1);
	  splice(@$refipot, $i, 1);
	  $i=0;
	}
      }
    }
    if ($$refkey[3] =~ /Y/){
      for ($i=0; $i<@$refx; $i++){
	if (abs($$refy[$i]) > $$refkey[4]){
	  splice(@$refx, $i, 1);
	  splice(@$refy, $i, 1);
	  splice(@$refz, $i, 1);
	  splice(@$refipot, $i, 1);
	  $i=0;
	}
      }
    }
    if ($$refkey[3] =~ /Z/){
      for ($i=0; $i<@$refx; $i++){
	if (abs($$refz[$i]) > $$refkey[4]){
	  splice(@$refx, $i, 1);
	  splice(@$refy, $i, 1);
	  splice(@$refz, $i, 1);
	  splice(@$refipot, $i, 1);
	  $i=0;
	}
      }
    }
  }
  if ($$refkey[0] =~ /Nx/){
    print "removing Neg X atoms\n";
    for ($i=0; $i<@$refx; $i++){
      if (($$refx[$i] < 0)){
	splice(@$refx, $i, 1);
	splice(@$refy, $i, 1);
	splice(@$refz, $i, 1);
	splice(@$refipot, $i, 1);
	$i=0;
      }
    }
  }
  if ($$refkey[0] =~ /Ny/){
    print "removing Neg Y atoms\n";
    for ($i=0; $i<@$refy; $i++){
      if (($$refy[$i] < 0)){
	splice(@$refx, $i, 1);
	splice(@$refy, $i, 1);
	splice(@$refz, $i, 1);
	splice(@$refipot, $i, 1);
	$i=0;
      }
    }
  }
  if ($$refkey[0] =~ /Nz/){
    print "removing Neg Z atoms\n";
    for ($i=0; $i<@$refz; $i++){
      if (($$refz[$i] < 0)){
	splice(@$refx, $i, 1);
	splice(@$refy, $i, 1);
	splice(@$refz, $i, 1);
	splice(@$refipot, $i, 1);
	$i=0;
      }
    }
  }
  $length = @$refx;
  print "number of atoms in cluster $length\n";
  return(@$refx,@$refy,@$refz,@$refipot);
}
